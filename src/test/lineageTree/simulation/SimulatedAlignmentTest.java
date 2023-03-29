package test.lineageTree.simulation;

import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.datatype.UserDataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.TreeParser;
import feast.simulation.SimulatedAlignment;
import lineageTree.simulation.SimulatedTypeWriterAlignment;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;
import org.junit.Before;
import org.junit.Test;

public class SimulatedAlignmentTest {

    // set up
    @Test
    public void testSimulationOnCherry(){

        Integer sequenceLength = 1;
        String outputFileName = "./simAl.nexus";
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        TypewriterSubstitutionModel submodel = new TypewriterSubstitutionModel();
        RealParameter insertrates = new RealParameter("5 5");
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("rates", insertrates, "frequencies", frequencies);
        submodel.calculateIntermediates();

        //site model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10");
        siteM.initByName( "gammaCategoryCount", 0,
                "substModel", submodel, "mutationRate", mutationRate);

        DataType integerData = new IntegerData();

        // simulate
        Randomizer.setSeed(1);
        SimulatedAlignment simAlignment = new SimulatedAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "sequenceLength", sequenceLength,
                "nrOfInsertionsPerSite", 100,
                "outputFileName", outputFileName,
                "userDataType", integerData
                );


    }

}
