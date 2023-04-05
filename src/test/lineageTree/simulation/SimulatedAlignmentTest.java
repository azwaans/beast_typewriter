package test.lineageTree.simulation;

import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import beast.util.TreeParser;
import feast.fileio.AlignmentFromNexus;
import lineageTree.simulation.SimulatedTypeWriterAlignment;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class SimulatedAlignmentTest {

    // set up
    @Test
    public void testSimulationOnCherry(){

        Integer sequenceLength = 1;
        String outputFileName = "./src/test/lineageTree/simulation/simAl.nexus";
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        TypewriterSubstitutionModel submodel = new TypewriterSubstitutionModel();
        RealParameter insertrates = new RealParameter("0.8 0.2");
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("editProbabilities", insertrates, "frequencies", frequencies);
        //submodel.setTargetBClength(1);


        //site model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10");
        siteM.initByName( "gammaCategoryCount", 0,
                "substModel", submodel, "mutationRate", mutationRate);

        DataType integerData = new IntegerData();

        // simulate
        Randomizer.setSeed(1);
        SimulatedTypeWriterAlignment simAlignment = new SimulatedTypeWriterAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "sequenceLength", sequenceLength,
                "nrOfInsertionsPerTarget", 100,
                "numberOfTargets", 1,
                "outputFileName", outputFileName,
                "userDataType", integerData
                );

        AlignmentFromNexus expectedAlignment = new AlignmentFromNexus();
        expectedAlignment.initByName("fileName",
                "src/test/lineageTree/simulation/expectedAlignment.nexus",
                "userDataType", integerData);


        assertEquals(expectedAlignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(expectedAlignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }

}
