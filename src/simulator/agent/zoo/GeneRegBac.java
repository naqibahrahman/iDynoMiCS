package simulator.agent.zoo;

import Jama.Matrix;
import odeSolver.GeneRegSolver;
import simulator.Simulator;
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
		out._proteinLevels = this._proteinLevels;
		return out;
	}
	
	public void initFromProtocolFile(Simulator aSim, XMLParser aSpeciesRoot)
	{
		super.initFromProtocolFile(aSim, aSpeciesRoot);
		
		_regulationSolver.init(_numProtTypes, getSpeciesParam().hmax, 
												getSpeciesParam().rtol);
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