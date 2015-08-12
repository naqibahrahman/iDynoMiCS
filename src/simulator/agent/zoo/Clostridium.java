package simulator.agent.zoo;

import java.util.ArrayList;

import Jama.Matrix;
import idyno.SimTimer;
import simulator.Simulator;
import simulator.SoluteGrid;
import simulator.agent.LocatedAgent;
import utils.ExtraMath;
import utils.XMLParser;

public class Clostridium extends GeneRegBac {

	private SoluteGrid _aipGrid;
	private SoluteGrid _acidGrid;
	private SoluteGrid _solventGrid;
	
	private Double _aipConc;
	private Double _acidConc;
	private Double _solventConc;
	
	private String _sporeStatus = "noSpore";
	
	public Clostridium()
	{
		super();
		
		_speciesParam = new ClostridiumParam();
		
		this._numProtTypes = 16;
		
		_proteinNames = new String[_numProtTypes];
		_proteinLevels = new Double[_numProtTypes];
		
		_proteinNames[0] = "SA";
		_proteinLevels[0] = 17500.0;
		
		_proteinNames[1] = "SAP";
		_proteinLevels[1] = 5000.0;
		
		_proteinNames[2] = "K";
		_proteinLevels[2] = 10.0;
		
		_proteinNames[3] = "KP";
		_proteinLevels[3] = 40.0;
		
		_proteinNames[4] = "Ph";
		_proteinLevels[4] = 50.0;
		
		_proteinNames[5] = "PhP";
		_proteinLevels[5] = 125.0;
		
		_proteinNames[6] = "Ab";
		_proteinLevels[6] = 6.0;
		
		_proteinNames[7] = "SigmaH";
		_proteinLevels[7] = 100.0;
		
		_proteinNames[8] = "A";
		_proteinLevels[8] = 1700.0;
		
		_proteinNames[9] = "AP";
		_proteinLevels[9] = 160.0;
		
		_proteinNames[10] = "B";
		_proteinLevels[10] = 40.0;
		
		_proteinNames[11] = "C";
		_proteinLevels[11] = 4.0;
		
		_proteinNames[12] = "S";
		_proteinLevels[12] = 0.01;
		
		_proteinNames[13] = "T";
		_proteinLevels[13] = 18500.0;
		
		_proteinNames[14] = "R";
		_proteinLevels[14] = 2.0;
		
		_proteinNames[15] = "RP";
		_proteinLevels[15] = 1850.0;
		
		
	}
	
	
	public Object clone() throws CloneNotSupportedException
	{
		Clostridium out = (Clostridium) super.clone();
		out._aipGrid = this._aipGrid;
		out._acidGrid = this._acidGrid;
		out._solventGrid = this._solventGrid;
		out._sporeStatus = this._sporeStatus;
		
		return out;
	}
	
	/**
	 * 
	 */
	@Override
	public void makeKid() throws CloneNotSupportedException
	{
		/*
		 * Create the new instance.
		 */
		Clostridium baby = (Clostridium) sendNewAgent();
		/*
		 * These are all generated randomly.
		 */
		this._myDivRadius = getDivRadius();
		baby._myDivRadius = getDivRadius();
		baby._myDeathRadius = getDeathRadius();
		/*
		 * Update the lineage.
		 */
		recordGenealogy(baby);
		/*
		 * Share mass of all compounds between two daughter cells and compute
		 * new size.
		 */
		divideCompounds(baby, getBabyMassFrac());
		/*
		 * In a chemostat, the daughter cells remain with the coordinates of
		 * their progenitor. Otherwise, compute movement to apply to both
		 * cells and apply it.
		 */
		if ( ! Simulator.isChemostat )
		{
			setDivisionDirection(getInteractDistance(baby)/2);
			baby._movement.subtract(_divisionDirection);
			_movement.add(_divisionDirection);
		}
		/*
		 * Now register the agent inside the guilds and the agent grid.
		 */
		baby.registerBirth();
		baby._netVolumeRate = 0.0;
		/*
		 * This is specific to Clostridium!
		 */
		if ( ExtraMath.getUniRandDbl() < getSpeciesParam().ProbParam)
			baby._sporeStatus = "canSpore";
		else
			baby._sporeStatus = "noSpore";
	}
	
	@Override
	protected void internalStep()
	{
		/*
		 * If it is a spore then we just return and do nothing
		 */
		if ( this._sporeStatus.equals("spore") )
			return;
		/*
		 * If cell is sporulating and the biomass of spore is bigger than the threshold
		 * then we turn off all reactions and set spore status as "spore"
		 */
		if ( this._sporeStatus.equals("sporulating") )
		{
			if ( this.particleMass[1] >= getSpeciesParam().SporeThresh) 
			{
				// turn off all reactions
		
				for ( int aReac : getSpeciesParam().offAllReactions )
					switchOffreaction(allReactions[aReac]);
				/*
				 * 
				 */
				this._sporeStatus = "spore";

				// set all protein to zero
				//for ( int i = 0; i < _proteinLevels.length; i++ )
				//	_proteinLevels[i] = 0.0;
				
						
				//turn on lysis
				
				
				return;
			}
			
		}
		/*
		 * If it is not a spore then we just continue doing other processes
		 */

		updateExternal();
		this._regulationSolver.setReferenceAgent(this);
		this._proteinLevels = this._regulationSolver.solve(this._proteinLevels, _agentGrid.AGENTTIMESTEP);
		checkSpo0A();
		secretesAIP();
		solventogenesis();
		updateExternal();
		/*
		 * Compute mass growth over all compartments.
		 */
		grow();
		/*
		 * Apply this mass growth of all compounds on global radius and mass.
		 */
		updateSize();
		/*
		 * Divide if you have to.
		 */
		if ( willDivide())
			divide();
		/*
		 * Die if you have to.
		 */
		if ( willDie() )
		{
			this.death = "tooSmall";
			die(true);
		}

	}
	
	
	@Override
	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot)
	{
		super.initFromProtocolFile(aSim, aSpeciesRoot);
		
		this._aipGrid = aSim.getSolute("aip");
		this._acidGrid = aSim.getSolute("acid");
		this._solventGrid = aSim.getSolute("solvent");
		
		
	}
	
	
	private void updateExternal()
	{
		_aipConc = this._aipGrid.getValueAt(this._location);
		_acidConc = this._acidGrid.getValueAt(this._location);
		_solventConc = this._solventGrid.getValueAt(this._location);
	}
	
	private Double activate(Double c, Double B, Double U, Double factor)
	{
		return c*B*factor / (( B * factor ) + U);	
	}
	
	private Double repress(Double c, Double B, Double U, Double factor)
	{
		return c*U / (( B * factor) + U);
	}
	
	private Double activateJac(Double c, Double B, Double U, Double factor)
	{
		return c*B*U / Math.pow(((B*factor)+U), 2);
	}
	private Double repressJac(Double c, Double B, Double U, Double factor)
	{
		return - c*B*U / Math.pow(((B*factor)+U), 2);
	}
	
	/**
	 * ODEs
	 *
	 */
	
	// ODE of dSA/dt
	private Double SA_Rate(Double SA, Double SAP, Double K, Double KP, Double Ph, Double SigmaH, Double AP)
	{
		Double rate = 0.0;
		rate = repress(getSpeciesParam().c_SA_l, getSpeciesParam().B_SAP_SA, getSpeciesParam().U_SAP_SA, SAP);
		rate += activate(1.0, getSpeciesParam().B_sigmaH_SA, getSpeciesParam().U_sigmaH_SA, SigmaH) * 
				(repress(getSpeciesParam().c_SA_2l, getSpeciesParam().B_SAP_SA, getSpeciesParam().U_SAP_SA, SAP) +
				activate(getSpeciesParam().c_SA_2h, getSpeciesParam().B_SAP_SA, getSpeciesParam().U_SAP_SA, SAP));
		rate -= 2 * getSpeciesParam().phi_AP_SA * AP * SA; 
		rate -= 2 * getSpeciesParam().phi_KP_SA * KP * SA; 
		rate += 2 * getSpeciesParam().phi_SAP_K * SAP * K; 
		rate += 2 * getSpeciesParam().psi_SAP * SAP; 
		rate += 2 * getSpeciesParam().phi_SAP_Ph * Ph * SAP; 
		rate -= getSpeciesParam().lambda_SA * SA; 
			
		return rate;
	}
		
	
	// ODE of dSAP/dt
	private Double SAP_Rate(Double SA, Double SAP, Double K, Double KP, Double Ph, Double AP)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().phi_KP_SA*KP*SA;
		rate -= getSpeciesParam().phi_SAP_K*SAP*K;
		rate += getSpeciesParam().phi_AP_SA*AP*SA;
		rate -= getSpeciesParam().phi_SAP_Ph*Ph*SAP;
		rate -= getSpeciesParam().psi_SAP*SAP;
		rate -= getSpeciesParam().lambda_SAP*SAP;
		return rate;
	}
		
	
	// ODE of dK/dt
	private Double K_Rate(Double SA, Double SAP, Double K, Double KP)
	{
		Double rate = 0.0;
		Double top = 0.0;
		Double bottom = 0.0;
			
		top = getSpeciesParam().c_K*_acidConc;
		bottom = _acidConc + 1;
		rate = top/bottom;
			
		rate -= getSpeciesParam().alpha*K;
		rate += getSpeciesParam().phi_KP_SA*KP*SA;
		rate -= getSpeciesParam().phi_SAP_K*SAP*K;
		rate += getSpeciesParam().psi_KP*KP;
		rate -= getSpeciesParam().lambda_K*K;
		return rate;
	}
	
	// ODE of dKP/dt
	private Double KP_Rate(Double SA, Double SAP, Double K, Double KP)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().alpha*K;
		rate -= getSpeciesParam().phi_KP_SA*KP*SA;
		rate += getSpeciesParam().phi_SAP_K*SAP*K;
		rate -= getSpeciesParam().psi_KP*KP;
		rate -= getSpeciesParam().lambda_KP*KP;
		return rate;
	}
	
	// ODE of dPh/dt
	private Double Ph_Rate(Double SAP, Double Ph, Double PhP)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().c_Ph;
		rate -= getSpeciesParam().phi_SAP_Ph*Ph*SAP;
		rate += getSpeciesParam().psi_PhP*PhP;
		rate -= getSpeciesParam().lambda_Ph*Ph;
		return rate;
	}

	// ODE of dPhP/dt
	private Double PhP_Rate(Double SAP, Double Ph, Double PhP)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().phi_SAP_Ph*Ph*SAP;
		rate -= getSpeciesParam().psi_PhP*PhP;
		rate -= getSpeciesParam().lambda_PhP*PhP;
		return rate;
	}
	
	// ODE of dAb/dt
	private Double Ab_Rate(Double SAP, Double Ab)
	{
		Double rate = 0.0;
		rate = 0.25*(repress(getSpeciesParam().c_Ab_1, getSpeciesParam().B_Ab_Ab, getSpeciesParam().U_Ab_Ab, Ab));
		rate += 0.25*(repress(getSpeciesParam().c_Ab_2, getSpeciesParam().B_SAP_Ab, getSpeciesParam().U_SAP_Ab, SAP)); 		
		rate -= getSpeciesParam().lambda_Ab*Ab;
		return rate;
	}
	
	// ODE of dSigmaH/dt
	private Double SigmaH_Rate(Double Ab, Double SigmaH)
	{
		Double rate = 0.0;
		rate = repress(getSpeciesParam().c_sigmaH, getSpeciesParam().B_Ab_sigmaH, getSpeciesParam().U_Ab_sigmaH, Ab );
		rate -= getSpeciesParam().lambda_sigmaH*SigmaH;
		return rate;
	}
	
	// ODE of dA/dt;
	private Double A_Rate(Double SA, Double A, Double AP, Double RP)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().c_agr_CA;
		rate -= getSpeciesParam().phi_RP_A*A*RP;
		rate += getSpeciesParam().phi_AP_SA*AP*SA;
		rate += getSpeciesParam().psi_AP*AP;
		rate -= getSpeciesParam().lambda_A*A;
		return rate;
	}
	
	// ODE of dAP/dt
	private Double AP_Rate(Double SA, Double A, Double AP, Double RP)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().phi_RP_A*A*RP;
		rate -= getSpeciesParam().phi_AP_SA*AP*SA;
		rate -= getSpeciesParam().psi_AP*AP;
		rate -= getSpeciesParam().lambda_AP*AP;
		return rate;
	}	
	
	// ODE of dB/dt
	private Double B_Rate(Double AP, Double B)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().c_agr_l;
		rate += activate(getSpeciesParam().c_agr_h, getSpeciesParam().B_AP_agr, getSpeciesParam().U_AP_agr, AP);
		rate -= getSpeciesParam().mu_agr*B;
		rate -= getSpeciesParam().lambda_B*B;
		return rate;	
	}
	
	//ODE of dC/dt
	private Double C_Rate(Double C)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().c_agr_CA;
		rate -= getSpeciesParam().mu_agr*C;
		rate -= getSpeciesParam().lambda_C*C;
		return rate;
			
	}
	
	//ODE of dS/dt
	private Double S_Rate(Double B, Double S, Double T)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().mu_agr*B;
		rate -= getSpeciesParam().k_agr*T*S;
		rate -= getSpeciesParam().lambda_S*S;
		return rate;
	}
	
	// ODE of dT/dt
	private Double T_Rate(Double B, Double T)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().mu_agr*B;
		rate -= getSpeciesParam().lambda_T*T;
		return rate;
	}
	
	// ODE of dR/dt
	private Double R_rate(Double C, Double R, Double RP)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().mu_agr*C;
		rate -= getSpeciesParam().beta_RP*R*_aipConc;
		rate += getSpeciesParam().gamma_RP*RP;
		rate -= getSpeciesParam().lambda_R*R;
		return rate;
	}
	
	// ODE of dRP/dt
	private Double RP_rate(Double R, Double RP)
	{
		Double rate = 0.0;
		rate = getSpeciesParam().beta_RP*R*_aipConc;
		rate -= getSpeciesParam().gamma_RP*RP;
		rate -= getSpeciesParam().lambda_RP*RP;
		return rate;
	}
	
	
	
	@Override
	public Matrix calc1stDeriv(Matrix levels)
	{
		Matrix rates = new Matrix(levels.getRowDimension(), 1);
		Double SA = levels.get(0, 0);
		Double SAP = levels.get(1, 0);
		Double K = levels.get(2, 0);
		Double KP = levels.get(3, 0);
		Double Ph = levels.get(4, 0);
		Double PhP = levels.get(5, 0);
		Double Ab = levels.get(6, 0);
		Double SigmaH = levels.get(7, 0);
		Double A = levels.get(8, 0);
		Double AP = levels.get(9, 0);
		Double B = levels.get(10, 0);
		Double C = levels.get(11, 0);
		Double S = levels.get(12, 0);
		Double T = levels.get(13, 0);
		Double R = levels.get(14, 0);
		Double RP = levels.get(15, 0);
		
		rates.set(0, 0, this.SA_Rate(SA, SAP, K, KP, Ph, SigmaH, AP));
		rates.set(1, 0, this.SAP_Rate(SA, SAP, K, KP, Ph, AP));
		rates.set(2, 0, this.K_Rate(SA, SAP, K, KP));
		rates.set(3, 0, this.KP_Rate(SA, SAP, K, KP));
		rates.set(4, 0, this.Ph_Rate(SAP, Ph, PhP));
		rates.set(5, 0, this.PhP_Rate(SAP, Ph, PhP));
		rates.set(6, 0, this.Ab_Rate(SAP, Ab));
		rates.set(7, 0, this.SigmaH_Rate(Ab, SigmaH));
		rates.set(8, 0, this.A_Rate(SAP, A, AP, RP));
		rates.set(9, 0, this.AP_Rate(SAP, A, AP, RP));
		rates.set(10, 0, this.B_Rate(AP, B));
		rates.set(11, 0, this.C_Rate(C));
		rates.set(12, 0, this.S_Rate(B, S, T));
		rates.set(13, 0, this.T_Rate(B, T));
		rates.set(14, 0, this.R_rate(C, R, RP));
		rates.set(15, 0, this.RP_rate(R, RP));
		
		
		return rates;
	}
	
	@Override
	public Matrix calcJacobian(Matrix levels)
	{
		Matrix dFdY = new Matrix(_numProtTypes, _numProtTypes, 0.0);
		
		Double SA = levels.get(0, 0);
		Double SAP = levels.get(1, 0);
		Double K = levels.get(2, 0);
		Double KP = levels.get(3, 0);
		Double Ph = levels.get(4, 0);
		Double Ab = levels.get(6, 0);
		Double SigmaH = levels.get(7, 0);
		Double A = levels.get(8, 0);
		Double AP = levels.get(9, 0);
		Double S = levels.get(12, 0);
		Double T = levels.get(13, 0);
		Double RP = levels.get(15, 0);
		
	/******************************************************************************************************	
		 * The following differentiate SA with respect to different proteins
	/******************************************************************************************************/
		
		//wrt SA
		dFdY.set(0, 0, -2 * getSpeciesParam().phi_AP_SA * AP - 2* getSpeciesParam().phi_KP_SA * KP - getSpeciesParam().lambda_SA);
		
		//wrt SAP
		dFdY.set(0, 1, repressJac(getSpeciesParam().c_SA_l, getSpeciesParam().B_SAP_SA, getSpeciesParam().U_SAP_SA, SAP) 
				+ (repressJac(getSpeciesParam().c_SA_2l, getSpeciesParam().B_SAP_SA, getSpeciesParam().U_SAP_SA, SAP) 
						* activate(1.0, getSpeciesParam().B_sigmaH_SA, getSpeciesParam().U_sigmaH_SA, SigmaH)) 
				+ (activateJac(getSpeciesParam().c_SA_2h, getSpeciesParam().B_SAP_SA, getSpeciesParam().U_SAP_SA, SAP) 
						* activate(1.0, getSpeciesParam().B_sigmaH_SA, getSpeciesParam().U_sigmaH_SA, SigmaH)) 
				+ 2* getSpeciesParam().phi_SAP_K*K + 2 * getSpeciesParam().psi_SAP + 2 * getSpeciesParam().phi_SAP_Ph*Ph) ;
		//wrt K
		dFdY.set(0, 2, 2*getSpeciesParam().phi_SAP_K*SAP);
	
		//wrt KP
		dFdY.set(0, 3, -2*getSpeciesParam().phi_KP_SA*SA);
		
		//wrt Ph
		dFdY.set(0, 4, 2*getSpeciesParam().phi_SAP_Ph*SAP);
		
		//wrt SigmaH
		dFdY.set(0, 7, activateJac(1.0, getSpeciesParam().B_sigmaH_SA, getSpeciesParam().U_sigmaH_SA, SigmaH) 
				* (repress(getSpeciesParam().c_SA_2l, getSpeciesParam().B_SAP_SA, getSpeciesParam().U_SAP_SA, SAP)
						+ activate(getSpeciesParam().c_SA_2h, getSpeciesParam().B_SAP_SA, getSpeciesParam().U_SAP_SA, SAP) ));
		//wrt AP
		dFdY.set(0, 9, -2*getSpeciesParam().phi_AP_SA*SA);
		
		
	/******************************************************************************************************	
		 * The following differentiate SAP with respect to different proteins
	/******************************************************************************************************/
		
		//wrt SA
		dFdY.set(1, 0,  getSpeciesParam().phi_KP_SA*KP + getSpeciesParam().phi_AP_SA*AP);
		
		//wrt SAP
		dFdY.set(1, 1, -getSpeciesParam().phi_SAP_K*K - getSpeciesParam().phi_SAP_Ph*Ph - getSpeciesParam().psi_SAP 
				+ getSpeciesParam().lambda_SAP);
	
		//wrt K
		dFdY.set(1, 2, -getSpeciesParam().phi_SAP_K*SAP);
		
		//wrt KP
		dFdY.set(1, 3, getSpeciesParam().phi_KP_SA*SA);
	
		//wrt Ph
		dFdY.set(1, 4, -getSpeciesParam().phi_SAP_Ph*SAP);
	
		//wrt AP
		dFdY.set(1, 9, getSpeciesParam().phi_AP_SA*SA);
	
		
	/******************************************************************************************************	
		 * The following differentiate K with respect to different proteins
	/******************************************************************************************************/
		
		//wrt SA
		dFdY.set(2, 0, getSpeciesParam().phi_KP_SA*KP);
	
		//wrt SAP
		dFdY.set(2, 1, -getSpeciesParam().phi_SAP_K*K);
		
		//wrt K
		dFdY.set(2, 2, -getSpeciesParam().phi_SAP_K*SAP + getSpeciesParam().lambda_K - getSpeciesParam().alpha);
		
		//wrt KP
		dFdY.set(2, 3, getSpeciesParam().phi_KP_SA*SA + getSpeciesParam().psi_KP);
		
		

	/******************************************************************************************************	
		 * The following differentiate KP with respect to different proteins
	/******************************************************************************************************/

		//wrt SA
		dFdY.set(3, 0, -getSpeciesParam().phi_KP_SA*KP);
			
		//wrt SAP
		dFdY.set(3, 1, getSpeciesParam().phi_SAP_K*K);
		
		//wrt K
		dFdY.set(3, 2, getSpeciesParam().alpha + getSpeciesParam().phi_SAP_K*SAP);
			
		//wrt KP
		dFdY.set(3, 3, -getSpeciesParam().phi_KP_SA*SA -getSpeciesParam().psi_KP -getSpeciesParam().lambda_KP);
	
		
	/******************************************************************************************************	
		 * The following differentiate Ph with respect to different proteins
	/******************************************************************************************************/
		
		//wrt SAP
		dFdY.set(4, 1, -getSpeciesParam().phi_SAP_Ph*Ph);
	
		//wrt Ph
		dFdY.set(4, 4, -getSpeciesParam().phi_SAP_Ph*SAP -getSpeciesParam().lambda_Ph);
	
		//wrt PhP
		dFdY.set(4, 5, getSpeciesParam().psi_PhP);
		
		
	/******************************************************************************************************	
		 * The following differentiate PhP with respect to different proteins
	/******************************************************************************************************/
		//wrt SAP
		dFdY.set(5, 1, getSpeciesParam().phi_SAP_Ph*Ph);
	
		//wrt Ph
		dFdY.set(5, 4, getSpeciesParam().phi_SAP_Ph*SAP);
		
		//wrt PhP
		dFdY.set(5, 5, -getSpeciesParam().psi_PhP -getSpeciesParam().lambda_PhP);
		
	/******************************************************************************************************	
		 * The following differentiate Ab with respect to different proteins
	/******************************************************************************************************/

		// wrt SAP
		dFdY.set(6, 1, 0.25 * (repressJac(getSpeciesParam().c_Ab_1, getSpeciesParam().B_Ab_Ab, getSpeciesParam().U_Ab_Ab, Ab) 
				+ repressJac(getSpeciesParam().c_Ab_2, getSpeciesParam().B_SAP_Ab, getSpeciesParam().U_SAP_Ab, SAP)));
		
		//wrt Ab
		dFdY.set(6, 6, -getSpeciesParam().lambda_Ab);
		
	/******************************************************************************************************	
		 * The following differentiate SigmaH with respect to different proteins
	/******************************************************************************************************/	
			
		//wrt Ab
		dFdY.set(7, 6, repressJac(getSpeciesParam().c_sigmaH, getSpeciesParam().B_Ab_sigmaH, getSpeciesParam().U_Ab_sigmaH, Ab));
	
		//wrt SigmaH
		dFdY.set(7, 7, -getSpeciesParam().lambda_sigmaH);
		
		
	/******************************************************************************************************	
		 * The following differentiate A with respect to different proteins
	/******************************************************************************************************/	
		
		//wrt SA
		dFdY.set(8, 0, getSpeciesParam().phi_AP_SA*AP);
			
		//wrt A
		dFdY.set(8, 8, -getSpeciesParam().phi_RP_A*RP -getSpeciesParam().lambda_A);
		
		//wrt AP
		dFdY.set(8, 9, getSpeciesParam().phi_AP_SA*SA + getSpeciesParam().psi_AP);
		
		//wrt RP
		dFdY.set(8, 15, -getSpeciesParam().phi_RP_A*A);

	/******************************************************************************************************	
		 * The following differentiate AP with respect to different proteins
	/******************************************************************************************************/	
		
		//wrt SA
		dFdY.set(9, 0, -getSpeciesParam().phi_AP_SA*AP);
	
		//wrt A
		dFdY.set(9, 8, getSpeciesParam().phi_RP_A*RP);
		
		//wrt AP
		dFdY.set(9, 9, -getSpeciesParam().phi_AP_SA*SA - getSpeciesParam().psi_AP - getSpeciesParam().lambda_AP);
	
		//wrt RP
		dFdY.set(9, 15, getSpeciesParam().phi_RP_A*A);
		
	/******************************************************************************************************	
		 * The following differentiate B with respect to different proteins
	/******************************************************************************************************/		
		
		//wrt AP
		dFdY.set(10, 9, activateJac(getSpeciesParam().c_agr_h, getSpeciesParam().B_AP_agr, getSpeciesParam().U_AP_agr, AP));
		
		//wrt B
		dFdY.set(10, 10, -getSpeciesParam().mu_agr - getSpeciesParam().lambda_B);
			
			
	/******************************************************************************************************	
		 * The following differentiate C with respect to different proteins
	/******************************************************************************************************/
			
		//wrt C
		dFdY.set(11, 11, -getSpeciesParam().mu_agr - getSpeciesParam().lambda_C);
		
	/******************************************************************************************************	
	 	 * The following differentiate S with respect to different proteins
	/******************************************************************************************************/
	
		//wrt B
		dFdY.set(12, 10, getSpeciesParam().mu_agr);
		
		//wrt S
		dFdY.set(12, 12, -getSpeciesParam().lambda_S - getSpeciesParam().k_agr*T);
		
		//wrt T
		dFdY.set(12, 13, -getSpeciesParam().k_agr*S);
		
		
	/******************************************************************************************************	
		 * The following differentiate T with respect to different proteins
	/******************************************************************************************************/
		//wrt B
		dFdY.set(13, 10, getSpeciesParam().mu_agr);
	
		//wrt T
		dFdY.set(13, 13, -getSpeciesParam().lambda_T);
		
	/******************************************************************************************************	
		 * The following differentiate R with respect to different proteins
	/******************************************************************************************************/
		
		//wrt C
		dFdY.set(14, 11, getSpeciesParam().mu_agr);
	
		//wrt R
		dFdY.set(14, 14, -getSpeciesParam().beta_RP*_aipConc - getSpeciesParam().lambda_R);
	
		//wrt RP
		dFdY.set(14, 15, getSpeciesParam().gamma_RP);
		
		
	/******************************************************************************************************	
		 * The following differentiate RP with respect to different proteins
	/******************************************************************************************************/
		
		//wrt R
		dFdY.set(15, 14, getSpeciesParam().beta_RP*_aipConc);
	
		//wrt RP
		dFdY.set(15, 15, -getSpeciesParam().gamma_RP - getSpeciesParam().lambda_RP);
	
		return dFdY;
	}
	
	public void secretesAIP()
	{
		Double S = _proteinLevels[12];
		Double T = _proteinLevels[13];
		Double R = _proteinLevels[14];
		Double RP = _proteinLevels[15];
		Double AIPrate = _aipConc;
		AIPrate += getSpeciesParam().k_agr*T*S;
		AIPrate -= getSpeciesParam().beta_RP*R;
		AIPrate += getSpeciesParam().gamma_RP*RP;
		
	// dilution 
		
		AIPrate *= this._volume/_aipGrid.getVoxelVolume();
		
		AIPrate *= SimTimer.getCurrentTimeStep();
		
		Double rate = Math.max(AIPrate, 0.0);
		
		_aipGrid.setValueAt(rate, this._location);
	}
	
	public void solventogenesis()
	{
		Double SAP = _proteinLevels[1]*this._volume/_aipGrid.getVoxelVolume();
		
		Double mu = getSpeciesParam().mu_i / (getSpeciesParam().mu_e + getSpeciesParam().mu_i);
		
		Double rate = mu * getSpeciesParam().rho * SAP * _acidConc;
		
		rate *= SimTimer.getCurrentTimeStep();
	
		Double solConc = Math.max(_solventConc + rate, 0.0);
		Double acConc = Math.max(_acidConc - rate, 0.0);
		
		_solventGrid.setValueAt(solConc, this._location);
		_acidGrid.setValueAt(acConc, this._location);
	}
	
	public void checkSpo0A()
	{
		Double SAP = this._proteinLevels[1];
		if ( SAP > getSpeciesParam().Spo0Athresh )
		{
			for ( int aReac : getSpeciesParam().offAllReactions )
				switchOffreaction(allReactions[aReac]);
			for ( int aReac : getSpeciesParam().onSolventogenesis )
				switchOnReaction(allReactions[aReac]);
			
			if ( this._sporeStatus.equals("canSpore") )
			{
				/*
				 * Turn on the reactions that should now be on.
				 */
				for ( int aReac : getSpeciesParam().onSporulation )
					switchOnReaction(allReactions[aReac]);
				
				
				this._sporeStatus = "sporulating";
			}
		}
	}


	public ClostridiumParam getSpeciesParam()
	{
		return (ClostridiumParam) _speciesParam;
	}
	
}