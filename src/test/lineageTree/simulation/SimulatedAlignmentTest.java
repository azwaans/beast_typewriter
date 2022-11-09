package test.lineageTree.simulation;

import beast.core.parameter.RealParameter;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.TreeParser;
import lineageTree.simulation.SimulatedAlignment;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;

public class SimulatedAlignmentTest {

    // set up
    public void testSimulationOnCherry(){

        Integer sequenceLength = 1;
        String outputFileName = "./simulatedAlignment.nexus";
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLAbelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        TypewriterSubstitutionModel submodel = new TypewriterSubstitutionModel();
        RealParameter insertrates = new RealParameter("0.5 0.5");
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("rates", insertrates, "frequencies", frequencies);
        submodel.calculateIntermediates();

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        // simulate
        Randomizer.setSeed(1);
        SimulatedAlignment simAlignment = new SimulatedAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "sequenceLength", sequenceLength,
                "outputFileName", outputFileName
                );


    }

}
