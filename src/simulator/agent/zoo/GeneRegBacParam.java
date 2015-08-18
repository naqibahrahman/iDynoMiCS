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
	 * 
	 */
	public Double rtol;
	
	public Double atol;
	
	/**
	 * 
	 */
	public Double hmax = 0.000001;
	
	
	
	
	public GeneRegBacParam()
	{
		super();
	}
	
	public void init(Simulator aSim, XMLParser aSpeciesRoot,
													XMLParser speciesDefaults)
	{
		super.init(aSim,aSpeciesRoot,speciesDefaults);
		
		rtol = getSpeciesParameterDouble("rtol", aSpeciesRoot, speciesDefaults);
		atol = getSpeciesParameterDouble("atol", aSpeciesRoot, speciesDefaults);
		if ( rtol == XMLParser.nullDbl && atol == XMLParser.nullDbl )
		{
			LogFile.writeLogAlways(
						"WARNING! No tolerance set in the chemostat solver.");
		}
		
		Double value = 
			getSpeciesParameterDouble("hmax", aSpeciesRoot, speciesDefaults);
		hmax = Double.isNaN(value) ? hmax : value;
		
	}
	
	
	
}