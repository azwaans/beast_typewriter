package sciphy;

import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.IntegerData;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import feast.fileio.AlignmentFromNexus;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

import sciphy.evolution.simulation.SimulatedSciPhyAlignment;
import sciphy.evolution.substitutionmodel.SciPhySubstitutionModel;

public class SimulatedAlignmentTest {

    // set up
    @Test
    public void testSimulationOnCherry() {

        Integer sequenceLength = 1;
        String outputFileName = "test/simAl.nexus";
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";

        Tree tree = new TreeParser();
        tree.initByName(
                "IsLabelledNewick", true,
                "newick", newick,
                "adjustTipHeights", false);

        SciPhySubstitutionModel submodel = new SciPhySubstitutionModel();
        RealParameter insertrates = new RealParameter("0.8 0.2");
        RealParameter freqs = new RealParameter("1.0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("editProbabilities", insertrates, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        RealParameter mutationRate = new RealParameter("10");
        siteM.initByName("gammaCategoryCount", 0,
                "substModel", submodel, "mutationRate", mutationRate);

        DataType integerData = new IntegerData();

        // simulate
        Randomizer.setSeed(1);
        SimulatedSciPhyAlignment simAlignment = new SimulatedSciPhyAlignment();
        simAlignment.initByName("tree", tree,
                "siteModel", siteM,
                "sequenceLength", sequenceLength,
                "arrayLength", 100,
                "numberOfTargets", 1,
                "outputFileName", outputFileName,
                "userDataType", integerData
        );

        AlignmentFromNexus expectedAlignment = new AlignmentFromNexus();
        expectedAlignment.initByName("fileName",
                "test/sciphy/expectedAlignment.nexus",
                "userDataType", integerData);


        assertEquals(expectedAlignment.getSequenceAsString("CHILD1"),
                simAlignment.getSequenceAsString("CHILD1"));

        assertEquals(expectedAlignment.getSequenceAsString("CHILD2"),
                simAlignment.getSequenceAsString("CHILD2"));

    }

}
