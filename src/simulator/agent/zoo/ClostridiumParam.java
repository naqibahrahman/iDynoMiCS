package simulator.agent.zoo;

import java.util.ArrayList;
import org.jdom.Element;

import simulator.Simulator;
import utils.XMLParser;

/**
 * 
 * 
 * @author Naqibah Rahman
 */

public class ClostridiumParam extends GeneRegBacParam
{
	/** 
	 * Clostridium default parameters
	 */
	/*
	 * Reaction indices
	 */
	protected ArrayList<Integer> onSporulation;
	protected ArrayList<Integer> offAllReactions;
	protected ArrayList<Integer> onSolventogenesis;
	
	
	/*
	 * Production 
	 */
	public Double c_H = 2.88e-8; //M cell-1 L h-1
	public Double c_K = 3.6e-07; //M h-1
	public Double c_agr_l = 1.44e-06; //M h-1
	public Double c_SA_l = 1.44e-06; //M h-1
	public Double c_SA_2l = 1.44e-06; //M h-1
	public Double c_Ph = 1.44e-06; //M h-1
	public Double c_Ab_1 = 1.44e-06; //M h-1
	public Double c_Ab_2 = 1.44e-06; //M h-1
	public Double c_sigmaH = 1.44e-06; //M h-1
	public Double k_agr = 1.44e+12; //M-1 h-1
	public Double c_agr_CA = 1.44e-05; //M h-1
	public Double c_agr_h = 0.000144; //M h-1
	public Double c_SA_2h = 0.000144; //M h-1
	/*
	 * DNA unbinding
	 */
	public Double U_SAP_SA = 2e-08; //M
	public Double U_sigmaH_SA = 2e-08; //M
	public Double U_AP_agr = 2e-08; //M
	public Double U_Ab_Ab = 2e-08; //M
	public Double U_SAP_Ab = 2e-08; //M
	public Double U_Ab_sigmaH = 2e-08; //M
	/*
	 * DNA binding
	 */
	public Double B_SAP_SA = 2e-08; //M
	public Double B_sigmaH_SA = 2e-08; //M
	public Double B_AP_agr = 2e-08; //M
	public Double B_Ab_Ab = 2e-08; //M
	public Double B_SAP_Ab = 2e-08; //M
	public Double B_Ab_sigmaH = 2e-08; //M
	/*
	 * Complex separation
	 */
	public Double gamma_RP = 360.0;  // h-1
	/*
	 * Complex formation
	 */
	public Double beta_RP = 298.8;  // h-1
	/*
	 * Degradation
	 */
	public Double lambda_SAP = 0.036;  // h-1
	public Double lambda_sigmaH = 1.44;  // h-1
	public Double lambda_A = 7.2;  // h-1
	public Double lambda_AP = 7.2;  // h-1
	public Double lambda_B = 7.2;  // h-1
	public Double lambda_S = 7.2;  // h-1
	public Double lambda_T = 7.2;  // h-1
	public Double lambda_R = 7.2;  // h-1
	public Double lambda_RP = 7.2;  // h-1
	public Double lambda_K = 7.2;  // h-1
	public Double lambda_KP = 7.2;  // h-1
	public Double lambda_SA = 7.2;  // h-1
	public Double lambda_Ab = 7.2;  // h-1
	public Double lambda_a = 72.0;  // h-1
	public Double lambda_Ph = 7.2;  // h-1
	public Double lambda_PhP = 7.2;  // h-1
	public Double lambda_C = 7.2;  // h-1
	
	/*
	 * Autophosphorylation
	 */
	public Double alpha = 360.0;  // h-1
	
	/*
	 * Phosphotransfer 
	 */
	public Double phi_KP_SA = 3600000.0; //M-1 h-1
	public Double phi_SAP_K = 3600000.0; //M-1 h-1
	public Double phi_AP_SA = 3600000.0; //M-1 h-1
	public Double phi_SAP_Ph = 3600000.0; //M-1 h-1
	public Double phi_RP_A = 3600000.0; //M-1 h-1
	
	/*
	 * Dephosphorylation 
	 */
	public Double psi_AP = 1.44; // h-1
	public Double psi_KP = 1.44; // h-1
	public Double psi_SAP = 1.44; // h-1
	public Double psi_PhP = 1.44; // h-1
	
	/*
	 * Transport
	 */
	public Double mu_e = 3600.0; //  h-1
	public Double mu_i = 3600.0; // h-1
	
	/*
	 * Movement to cell membrane
	 */
	public Double mu_agr = 3600.0; // h-1
	
	/**
	 * Threshold intracellular Spo0A~P concentration triggering
	 * solventogenesis and/or sporulation.
	 */
	public Double spo0APthresh = 10e-06; //M
	
	/**
	 * Threshold spore mass triggering the end of sporulation.
	 */
	public Double sporeThresh = 6000.0; // fg
	
	/**
	 * Probability that a newly-created daughter cell will sporulate if the 
	 * spo0A~P threshold is ever reached.
	 */
	public Double sporeProb = 0.5;
	
	public ClostridiumParam()
	{
		super();
	}
	
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		super.init(aSim,aSpeciesRoot,speciesDefaults);
		Double value;
		
		// Production
		
		value = getSpeciesParameterDouble("c_H", aSpeciesRoot, speciesDefaults);
		c_H = Double.isNaN(value) ? c_H : value;
		
		value = getSpeciesParameterDouble("c_K", aSpeciesRoot, speciesDefaults);
		c_K = Double.isNaN(value) ? c_K : value;
		
		value = getSpeciesParameterDouble("c_agr_l", aSpeciesRoot, speciesDefaults);
		c_agr_l = Double.isNaN(value) ? c_agr_l : value;
		
		value = getSpeciesParameterDouble("c_SA_l", aSpeciesRoot, speciesDefaults);
		c_SA_l = Double.isNaN(value) ? c_SA_l : value;
		
		value = getSpeciesParameterDouble("c_SA_2l", aSpeciesRoot, speciesDefaults);
		c_SA_2l = Double.isNaN(value) ? c_SA_2l : value;
		
		value = getSpeciesParameterDouble("c_Ph", aSpeciesRoot, speciesDefaults);
		c_Ph = Double.isNaN(value) ? c_Ph : value;
		
		value = getSpeciesParameterDouble("c_Ab_1", aSpeciesRoot, speciesDefaults);
		c_Ab_1 = Double.isNaN(value) ? c_Ab_1 : value;
		
		value = getSpeciesParameterDouble("c_Ab_2", aSpeciesRoot, speciesDefaults);
		c_Ab_2 = Double.isNaN(value) ? c_Ab_2 : value;
		
		value = getSpeciesParameterDouble("c_sigmaH", aSpeciesRoot, speciesDefaults);
		c_sigmaH = Double.isNaN(value) ? c_sigmaH : value;
		
		value = getSpeciesParameterDouble("k_agr", aSpeciesRoot, speciesDefaults);
		k_agr = Double.isNaN(value) ? k_agr : value;
		
		value = getSpeciesParameterDouble("c_agr_CA", aSpeciesRoot, speciesDefaults);
		c_agr_CA = Double.isNaN(value) ? c_agr_CA : value;
		
		value = getSpeciesParameterDouble("c_agr_h", aSpeciesRoot, speciesDefaults);
		c_agr_h = Double.isNaN(value) ? c_agr_h : value;
		
		value = getSpeciesParameterDouble("c_SA_2h", aSpeciesRoot, speciesDefaults);
		c_SA_2h = Double.isNaN(value) ? c_SA_2h : value;
		
		// DNA unbinding
		
		value = getSpeciesParameterDouble("U_SAP_SA", aSpeciesRoot, speciesDefaults);
		U_SAP_SA = Double.isNaN(value) ? U_SAP_SA : value;
		
		value = getSpeciesParameterDouble("U_sigmaH_SA", aSpeciesRoot, speciesDefaults);
		U_sigmaH_SA = Double.isNaN(value) ? U_sigmaH_SA : value;
		
		value = getSpeciesParameterDouble("U_AP_agr", aSpeciesRoot, speciesDefaults);
		U_AP_agr = Double.isNaN(value) ? U_AP_agr : value;
		
		value = getSpeciesParameterDouble("U_Ab_Ab", aSpeciesRoot, speciesDefaults);
		U_Ab_Ab = Double.isNaN(value) ? U_Ab_Ab : value;
		
		value = getSpeciesParameterDouble("U_SAP_Ab", aSpeciesRoot, speciesDefaults);
		U_SAP_Ab = Double.isNaN(value) ? U_SAP_Ab : value;
		
		value = getSpeciesParameterDouble("U_Ab_sigmaH", aSpeciesRoot, speciesDefaults);
		U_Ab_sigmaH = Double.isNaN(value) ? U_Ab_sigmaH : value;
		
		// DNA binding
		
		value = getSpeciesParameterDouble("B_SAP_SA", aSpeciesRoot, speciesDefaults);
		B_SAP_SA = Double.isNaN(value) ? B_SAP_SA : value;
		
		value = getSpeciesParameterDouble("B_sigmaH_SA", aSpeciesRoot, speciesDefaults);
		B_sigmaH_SA = Double.isNaN(value) ? B_sigmaH_SA : value;
		
		value = getSpeciesParameterDouble("B_AP_agr", aSpeciesRoot, speciesDefaults);
		B_AP_agr = Double.isNaN(value) ? B_AP_agr : value;
		
		value = getSpeciesParameterDouble("B_Ab_Ab", aSpeciesRoot, speciesDefaults);
		B_Ab_Ab = Double.isNaN(value) ? B_Ab_Ab : value;
		
		value = getSpeciesParameterDouble("B_SAP_Ab", aSpeciesRoot, speciesDefaults);
		B_SAP_Ab = Double.isNaN(value) ? B_SAP_Ab : value;
		
		value = getSpeciesParameterDouble("B_Ab_sigmaH", aSpeciesRoot, speciesDefaults);
		B_Ab_sigmaH = Double.isNaN(value) ? B_Ab_sigmaH : value;
		
		// Complex separation
		
		value = getSpeciesParameterDouble("gamma_RP", aSpeciesRoot, speciesDefaults);
		gamma_RP = Double.isNaN(value) ? gamma_RP : value;
		
		// Complex formation
		
		value = getSpeciesParameterDouble("beta_RP", aSpeciesRoot, speciesDefaults);
		beta_RP = Double.isNaN(value) ? beta_RP : value;
		
		//Degradation
		
		value = getSpeciesParameterDouble("lambda_SAP", aSpeciesRoot, speciesDefaults);
		lambda_SAP = Double.isNaN(value) ? lambda_SAP : value;
		
		value = getSpeciesParameterDouble("lambda_sigmaH", aSpeciesRoot, speciesDefaults);
		lambda_sigmaH = Double.isNaN(value) ? lambda_sigmaH : value;
		
		value = getSpeciesParameterDouble("lambda_A", aSpeciesRoot, speciesDefaults);
		lambda_A = Double.isNaN(value) ? lambda_A : value;
		
		value = getSpeciesParameterDouble("lambda_AP", aSpeciesRoot, speciesDefaults);
		lambda_AP = Double.isNaN(value) ? lambda_AP : value;
		
		value = getSpeciesParameterDouble("lambda_B", aSpeciesRoot, speciesDefaults);
		lambda_B = Double.isNaN(value) ? lambda_B : value;
		
		value = getSpeciesParameterDouble("lambda_S", aSpeciesRoot, speciesDefaults);
		lambda_S = Double.isNaN(value) ? lambda_S : value;
		
		value = getSpeciesParameterDouble("lambda_T", aSpeciesRoot, speciesDefaults);
		lambda_T = Double.isNaN(value) ? lambda_T : value;
		
		value = getSpeciesParameterDouble("lambda_R", aSpeciesRoot, speciesDefaults);
		lambda_R = Double.isNaN(value) ? lambda_R : value;
		
		value = getSpeciesParameterDouble("lambda_RP", aSpeciesRoot, speciesDefaults);
		lambda_RP = Double.isNaN(value) ? lambda_RP : value;
		
		value = getSpeciesParameterDouble("lambda_K", aSpeciesRoot, speciesDefaults);
		lambda_K = Double.isNaN(value) ? lambda_K : value;
		
		value = getSpeciesParameterDouble("lambda_KP", aSpeciesRoot, speciesDefaults);
		lambda_KP = Double.isNaN(value) ? lambda_KP : value;
		
		value = getSpeciesParameterDouble("lambda_SA", aSpeciesRoot, speciesDefaults);
		lambda_SA = Double.isNaN(value) ? lambda_SA : value;
		
		value = getSpeciesParameterDouble("lambda_Ab", aSpeciesRoot, speciesDefaults);
		lambda_Ab = Double.isNaN(value) ? lambda_Ab : value;
		
		value = getSpeciesParameterDouble("lambda_a", aSpeciesRoot, speciesDefaults);
		lambda_a = Double.isNaN(value) ? lambda_a : value;
		
		value = getSpeciesParameterDouble("lambda_Ph", aSpeciesRoot, speciesDefaults);
		lambda_Ph = Double.isNaN(value) ? lambda_Ph : value;
		
		
		value = getSpeciesParameterDouble("lambda_PhP", aSpeciesRoot, speciesDefaults);
		lambda_PhP = Double.isNaN(value) ? lambda_PhP : value;
		
		
		value = getSpeciesParameterDouble("lambda_C", aSpeciesRoot, speciesDefaults);
		lambda_C = Double.isNaN(value) ? lambda_C : value;
		
		// Autophosphorylation
		
		value = getSpeciesParameterDouble("alpha", aSpeciesRoot, speciesDefaults);
		alpha = Double.isNaN(value) ? alpha : value;
		
		// Phosphotransfer
		
		value = getSpeciesParameterDouble("phi_KP_SA", aSpeciesRoot, speciesDefaults);
		phi_KP_SA = Double.isNaN(value) ? phi_KP_SA : value;
		
		
		value = getSpeciesParameterDouble("phi_SAP_K", aSpeciesRoot, speciesDefaults);
		phi_SAP_K = Double.isNaN(value) ? phi_SAP_K : value;
		
		
		value = getSpeciesParameterDouble("phi_AP_SA", aSpeciesRoot, speciesDefaults);
		phi_AP_SA = Double.isNaN(value) ? phi_AP_SA : value;
		
		
		value = getSpeciesParameterDouble("phi_SAP_Ph", aSpeciesRoot, speciesDefaults);
		phi_SAP_Ph = Double.isNaN(value) ? phi_SAP_Ph : value;
		
		
		value = getSpeciesParameterDouble("phi_RP_A", aSpeciesRoot, speciesDefaults);
		phi_RP_A = Double.isNaN(value) ? phi_RP_A : value;
		
		// Dephosphorylation
		
		value = getSpeciesParameterDouble("psi_AP", aSpeciesRoot, speciesDefaults);
		psi_AP = Double.isNaN(value) ? psi_AP : value;
		
		
		value = getSpeciesParameterDouble("psi_KP", aSpeciesRoot, speciesDefaults);
		psi_KP = Double.isNaN(value) ? psi_KP : value;
		
		
		value = getSpeciesParameterDouble("psi_SAP", aSpeciesRoot, speciesDefaults);
		psi_SAP = Double.isNaN(value) ? psi_SAP : value;
		
		
		value = getSpeciesParameterDouble("psi_PhP", aSpeciesRoot, speciesDefaults);
		psi_PhP = Double.isNaN(value) ? psi_PhP : value;
		
		// Transport
		
		value = getSpeciesParameterDouble("mu_e", aSpeciesRoot, speciesDefaults);
		mu_e = Double.isNaN(value) ? mu_e : value;
		
		
		value = getSpeciesParameterDouble("mu_i", aSpeciesRoot, speciesDefaults);
		mu_i = Double.isNaN(value) ? mu_i : value;
		
		// Movement to cell membrane
		
		value = getSpeciesParameterDouble("mu_agr", aSpeciesRoot, speciesDefaults);
		mu_agr = Double.isNaN(value) ? mu_agr : value;
		
		
		// Acid conversion
		
		//value = getSpeciesParameterDouble("rho", aSpeciesRoot, speciesDefaults);
		//rho = Double.isNaN(value) ? rho : value;
		
		//Threshold
		value = getSpeciesParameterDouble("Spo0Athresh", aSpeciesRoot, speciesDefaults);
		spo0APthresh = Double.isNaN(value) ? spo0APthresh : value;
		
		value = getSpeciesParameterDouble("SporeThresh", aSpeciesRoot, speciesDefaults);
		sporeThresh = Double.isNaN(value) ? sporeThresh : value;
		
		//Probability
		value = getSpeciesParameterDouble("sporeProb", aSpeciesRoot, speciesDefaults);
		sporeProb = Double.isNaN(value) ? sporeProb : value;
		
		
		
		/*
		 * Sporulation switch
		 */
		int reacIndex;
		onSporulation = new ArrayList<Integer>();
		
		XMLParser switchParser = new XMLParser(aSpeciesRoot.getChildElement("reactionSwitch"));
		XMLParser childParser;
		
		//////////////////////////////////////////////////////////////////////////
		// create list of reactions during the On state (sporulation)
		childParser = new XMLParser(switchParser.getChildElement("sporulating"));
		for (Element aReactionMarkUp : childParser.getChildrenElements("reaction"))
		{
			// Add the reaction to the list of reactions for the on-state
			// but add only ACTIVE reactions!!
			reacIndex = aSim.getReactionIndex(aReactionMarkUp.getAttributeValue("name"));
			onSporulation.add(reacIndex);
		}
		
		//////////////////////////////////////////////////////////////////////////
		// create list of reactions during the Off state
		offAllReactions = new ArrayList<Integer>();
		for ( int iReac = 0; iReac < aSim.reactionList.length; iReac++)
			offAllReactions.add(iReac);
		
		//////////////////////////////////////////////////////////////////////////
		// create list of reactions during the On state (solventogenesis)
		onSolventogenesis = new ArrayList<Integer>();
		childParser = new XMLParser(switchParser.getChildElement("solventogenesis"));
		for (Element aReactionMarkUp : childParser.getChildrenElements("reaction"))
		{
		// Add the reaction to the list of reactions for the on-state
			// but add only ACTIVE reactions!!
			reacIndex = aSim.getReactionIndex(aReactionMarkUp.getAttributeValue("name"));
			onSolventogenesis.add(reacIndex);
		}
	
		
	}
}
