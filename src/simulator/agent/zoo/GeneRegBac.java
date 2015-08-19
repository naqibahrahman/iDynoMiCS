package simulator.agent.zoo;

import Jama.Matrix;
import odeSolver.GeneRegSolver;
import simulator.Simulator;
import simulator.geometry.ContinuousVector;
import utils.ExtraMath;
import utils.LogFile;
import utils.XMLParser;

public abstract class GeneRegBac extends Bacterium implements Cloneable
{
	/**
	 * Ordinary Differential Solver used to calculate the changes in gene
	 * expression.
	 */
	protected GeneRegSolver _regulationSolver = new GeneRegSolver();
	
	/**
	 * 
	 */
	protected int _numProtTypes;
	
	/**
	 * 
	 */
	protected String[] _proteinNames;
	
	/**
	 * 
	 */
	protected Double[] _proteinLevels;
	
	
	
	public GeneRegBac()
	{
		super();
		_speciesParam = new GeneRegBacParam();
	}
	
	
	public Object clone() throws CloneNotSupportedException
	{
		GeneRegBac out = (GeneRegBac) super.clone();
		out._regulationSolver = this._regulationSolver;
		out._numProtTypes = this._numProtTypes;
		out._proteinNames = this._proteinNames;
		out._proteinLevels = new Double[this._numProtTypes];
		for ( int i = 0; i < this._numProtTypes; i++ )
			out._proteinLevels[i] = this._proteinLevels[i];
		return out;
	}
	
	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot)
	{
		super.initFromProtocolFile(aSim, aSpeciesRoot);
		
		_regulationSolver.init(_numProtTypes, getSpeciesParam().hmax, 
							getSpeciesParam().rtol, getSpeciesParam().atol);
		_regulationSolver.setReferenceAgent(this);
	}
	
	@Override
	public void initFromResultFile(Simulator aSim, String[] singleAgentData) 
	{
		/*
		 * 
		 */
		int iDataStart = singleAgentData.length - this._numProtTypes;
		for ( int i = 0; i < this._numProtTypes; i++ )
			this._proteinLevels[i] = Double.parseDouble(singleAgentData[iDataStart+i]);
		/*
		 * Now go up the hierarchy with the rest of the data.
		 */
		String[] remainingSingleAgentData = new String[iDataStart];
		for (int i=0; i<iDataStart; i++)
			remainingSingleAgentData[i] = singleAgentData[i];
		super.initFromResultFile(aSim, remainingSingleAgentData);
	}
	
	/**
	 * \brief Create a new agent in a specified position.
	 * 
	 * <p>Copied from the same method in LocatedAgent, but with the method
	 * {@link #randomiseProteinLevels()} added.</p>
	 * 
	 * @param position	Vector stating where this agent should be located.
	 */
	@Override
	public void createNewAgent(ContinuousVector position) 
	{
		System.out.println("Progenitor "+_proteinNames[5]+" is "+_proteinLevels[5]); //BUGHUNT
		try 
		{
			// Get a clone of the progenitor.
			GeneRegBac baby = (GeneRegBac) sendNewAgent();
			baby.giveName();
			baby.updateSize();
			
			this._myDivRadius = getDivRadius();
			baby._myDivRadius = getDivRadius();
			baby._myDeathRadius = getDeathRadius();
			
			// Just to avoid to be in the carrier.
			// TODO Rob 13Mar2015: Is this correct?
			position.x += this._totalRadius;
			
			baby.setLocation(position);
			baby.registerBirth();
			baby.randomiseProteinLevels();
		} 
		catch (CloneNotSupportedException e) 
		{
			LogFile.writeError(e, "GeneRegBac.createNewAgent()");
		}
	}
	
	protected void randomiseProteinLevels()
	{
		double cv = getSpeciesParam().initialProteinCV;
		if ( cv == 0.0 )
			return;
		for ( int i = 0; i < _numProtTypes; i++ )
			_proteinLevels[i] = ExtraMath.deviateFromCV(_proteinLevels[i], cv);
	}
	
	public abstract Matrix calc1stDeriv(Matrix levels);
	
	public abstract Matrix calcJacobian(Matrix levels);
	
	
	public GeneRegBacParam getSpeciesParam()
	{
		return (GeneRegBacParam) _speciesParam;
	}
	
	
	@Override
	public StringBuffer sendHeader()
	{
		// return the header file for this agent's values after sending those for super
		StringBuffer tempString = super.sendHeader();
		
		//all proteinNames - use loop
		for ( int i = 0; i < this._numProtTypes; i++ )
			tempString.append("," + this._proteinNames[i]);
		return tempString;
	}
	
	@Override
	public StringBuffer writeOutput()
	{
		// write the data matching the header file
		StringBuffer tempString = super.writeOutput();
		
		for ( int i = 0; i < this._numProtTypes; i++ )
			tempString.append("," + this._proteinLevels[i]);
		return tempString;
	}
}