package test;

import beast.core.Description;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import lineageTree.distributions.TypewriterTreeLikelihood;
import lineageTree.substitutionmodel.TypewriterSubstitutionModelHomogeneous;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.junit.Before;
import org.junit.Test;

import java.util.List;


public class TypewriterLikelihoodTest {

    TypewriterTreeLikelihood typewriterLikelihood;
    TypewriterSubstitutionModelHomogeneous substModel;
    SiteModel siteM;
    StrictClockModel clockModel;
    RealParameter origin;
    IntegerParameter arraylength;

    @Before
    public void setUp(){
        typewriterLikelihood = new TypewriterTreeLikelihood();
        substModel = new TypewriterSubstitutionModelHomogeneous();
        TypewriterSubstitutionModelHomogeneous substitutionModel = new TypewriterSubstitutionModelHomogeneous();
        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        substitutionModel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);

        //site model

        siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", substitutionModel);


        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);

         origin = new RealParameter("6");
         arraylength = new IntegerParameter("5");

    }

    @Test (expected = RuntimeException.class)
    public void testExceptionForSingleBranchInput(){

        String newick = "(CHILD1:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");

        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        // this should fail
        typewriterLikelihood.initByName("data",alignment,"tree",tree,"siteModel",siteM,"branchRateModel",clockModel,"origin", origin,"arrayLength", arraylength);
    }

    @Test (expected = RuntimeException.class)
    public void testExceptionForSingleNodeInput(){

        String newick = "(CHILD1:0)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");

        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);


        // this should fail
        typewriterLikelihood.initByName("data",alignment,"tree",tree,"siteModel",siteM,"branchRateModel",clockModel,"origin", origin,"arrayLength", arraylength);
    }

    @Test (expected = RuntimeException.class)
    public void testExceptionForOriginNegativeInput(){
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");


        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);


        // this should fail
        RealParameter originNegative = new RealParameter("-1");
        typewriterLikelihood.initByName("data",alignment,"tree",tree,"siteModel",siteM,"branchRateModel",clockModel,"origin", originNegative,"arrayLength", arraylength);

    }

    @Test (expected = RuntimeException.class)
    public void testExceptionForOriginSmallerThanTreeHeightInput(){
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");


        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        RealParameter originNegative = new RealParameter("4");
        // this should fail
        typewriterLikelihood.initByName("data",alignment,"tree",tree,"siteModel",siteM,"branchRateModel",clockModel,"origin", originNegative,"arrayLength", arraylength);
        double LogPCalc = typewriterLikelihood.calculateLogP();


    }

    @Test (expected = RuntimeException.class)
    public void testExceptionForNegativeTargetLengthInput(){
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");


        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);


        // this should fail
        IntegerParameter arrayNegative = new IntegerParameter("-1");
        typewriterLikelihood.initByName("data",alignment,"tree",tree,"siteModel",siteM,"branchRateModel",clockModel,"origin", origin,"arrayLength", arrayNegative);

    }
    //TODO maybe automatically set the array length based on input data?
    @Test (expected = RuntimeException.class)
    public void testExceptionForMismatchedTargetLengthInput(){
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");


        Tree tree = new TreeParser();
        tree.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);


        // this should fail
        IntegerParameter arrayTooLarge = new IntegerParameter("6");

        typewriterLikelihood.initByName("data",alignment,"tree",tree,"siteModel",siteM,"branchRateModel",clockModel,"origin", origin,"arrayLength", arrayTooLarge);

    }
}
