package simulator.agent.zoo;

import simulator.Simulator;
import utils.LogFile;
import utils.XMLParser;

/**
 * 
 * 
 * @author Robert Clegg (rjc096@bham.ac.uk), Centre for Systems Biology,
 * University of Birmingham (UK).
 */
public class GeneRegBacParam extends BacteriumParam 
{
	

	/**
	 * Relative tolerance of the ODE solver. rtol and/or atol must be set.
	 */
	public Double rtol;
	
	/**
	 * Absolute tolerance of the ODE solver. rtol and/or atol must be set.
	 */
	public Double atol;
	
	/**
	 * Maximum step size of the ODE solver.
	 */
	public Double hmax = 1e-6;
	
	/**
	 * Coefficient of Variation of the protein levels for a cell created via
	 * the protocol file.  
	 */
	public Double initialProteinCV = 0.0;
	
	
	public GeneRegBacParam()
	{
		super();
	}
	
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		super.init(aSim,aSpeciesRoot,speciesDefaults);
		
		rtol = getSpeciesParameterDouble("rtol",aSpeciesRoot,speciesDefaults);
		atol = getSpeciesParameterDouble("atol",aSpeciesRoot,speciesDefaults);
		if ( rtol == XMLParser.nullDbl && atol == XMLParser.nullDbl )
		{
			LogFile.writeLogAlways(
						"WARNING! No tolerance set in the chemostat solver.");
		}
		
		Double value;
		
		value =	getSpeciesParameterDouble("hmax",
											aSpeciesRoot, speciesDefaults);
		hmax = Double.isNaN(value) ? hmax : value;
		
		value = getSpeciesParameterDouble("initialProteinCV",
											aSpeciesRoot, speciesDefaults);
		initialProteinCV = Double.isNaN(value) ? initialProteinCV : value;
		
	}
	
	
	
}