package test;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.*;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.datatype.DataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
//import beast.evolution.datatype.integer;
import beast.evolution.tree.Node;
import lineageTree.distributions.TypewriterTreeLikelihood;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import lineageTree.substitutionmodel.TypewriterSubstitutionModelHomogeneous;
import org.junit.Test;

import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.*;


@Description("Tests for Typewriter Data on the ancestral state reconstruction and the likelihood")
public class TypewriterDataTest {

@Test
public void test_typewriter_data() {
    Sequence a = new Sequence("cell1", "2,1,0,0,0");
    Sequence b = new Sequence("cell2", "1,2,0,0,0");

    Alignment alignment = new Alignment();
    alignment.initByName("sequence", a, "dataType", "integer");
    alignment.initByName("sequence", b, "dataType", "integer");

    //internal representation of the sequences for the package:
    List<Integer> sequence_a = alignment.getCounts().get(0);
    List<Integer> sequence_b = alignment.getCounts().get(1);

    //ancestral sequences
    List<List<Integer>> ancs_sequence_a = TypewriterTreeLikelihood.get_possible_ancestors(sequence_a);
    List<List<Integer>> ancs_sequence_b = TypewriterTreeLikelihood.get_possible_ancestors(sequence_b);

    //manually create ancestral states
    List<Integer> allele21 = Arrays.asList(2, 1, 0, 0, 0);
    List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
    List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
    List<Integer> allele2 = Arrays.asList(2, 0, 0, 0, 0);
    List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

    // For a sequence with n sites, there are n_edited_sites + 1 ancestral sequences (all edited position + fully unedited))
    assertEquals(ancs_sequence_a.size(), 3, 1e-5);
    Log.info.println("aancers" + ancs_sequence_a);
    //for any sequence, its possible ancestral sequences are itself + removing edits 1 by 1 + unedited

    //check for sequence a
    assertTrue(ancs_sequence_a.contains(allele21));
    assertTrue(ancs_sequence_a.contains(allele2));
    assertTrue(ancs_sequence_a.contains(allele0));

    //check for sequence b
    assertTrue(ancs_sequence_b.contains(allele12));
    assertTrue(ancs_sequence_b.contains(allele1));
    assertTrue(ancs_sequence_b.contains(allele0));



}

    @Test
    public void test_ancestors_identical_sequences() {
        // Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();
        RealParameter freqs = new RealParameter("1.0 0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel);

        //test ancestral states sets calculations
        likelihood.traverseAncestral(tree1.getRoot());
        Hashtable<Integer,List<List<Integer>>> statesDictionary = likelihood.ancestralStates;


        //first calculate states dictionary
        //manually create states:
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(0).size(),3);
        assertTrue(statesDictionary.get(0).contains(allele12));
        assertTrue(statesDictionary.get(0).contains(allele1));
        assertTrue(statesDictionary.get(0).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(1).size(),3);
        assertTrue(statesDictionary.get(1).contains(allele12));
        assertTrue(statesDictionary.get(1).contains(allele1));
        assertTrue(statesDictionary.get(1).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(2).size(),3);
        assertTrue(statesDictionary.get(2).contains(allele12));
        assertTrue(statesDictionary.get(2).contains(allele1));
        assertTrue(statesDictionary.get(2).contains(allele0));



    }

    @Test
    public void test_ancestors() {
        // Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,3,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter freqs = new RealParameter("1.0 0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName( "editfrequencies", freqs, "frequencies", frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel);

        //test ancestral states sets calculations
        likelihood.traverseAncestral(tree1.getRoot());
        Hashtable<Integer,List<List<Integer>>> statesDictionary = likelihood.ancestralStates;


        //first calculate states dictionary
        //manually create states:
        List<Integer> allele123 = Arrays.asList(1, 2, 3, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(3, statesDictionary.size());

        //1st leaf
        assertEquals(statesDictionary.get(0).size(),4);
        assertTrue(statesDictionary.get(0).contains(allele123));
        assertTrue(statesDictionary.get(0).contains(allele12));
        assertTrue(statesDictionary.get(0).contains(allele1));
        assertTrue(statesDictionary.get(0).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(1).size(),3);
        assertTrue(statesDictionary.get(1).contains(allele12));
        assertTrue(statesDictionary.get(1).contains(allele1));
        assertTrue(statesDictionary.get(1).contains(allele0));

        //root node
        assertEquals(statesDictionary.get(2).size(),3);
        assertTrue(statesDictionary.get(2).contains(allele12));
        assertTrue(statesDictionary.get(2).contains(allele1));
        assertTrue(statesDictionary.get(2).contains(allele0));



    }




    @Test
    public void test_ancestral_sets_3leaf() {


        // Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "((A:1,B:1):1,C:2)";
        Sequence a = new Sequence("A", "1,2,0,0,0");
        Sequence b = new Sequence("B", "1,1,0,0,0");
        Sequence c = new Sequence("C", "2,1,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");
        alignment.initByName("sequence", c, "dataType", "integer");


        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter freqs = new RealParameter("0.5 0.5");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("editfrequencies", freqs,"frequencies",frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel);


        //calculate states dictionary
        likelihood.traverseAncestral(tree1.getRoot());
        Hashtable<Integer,List<List<Integer>>> statesDictionary = likelihood.ancestralStates;
        assertEquals(5, statesDictionary.size());

        //Manually create states
        List<Integer> allele12 = Arrays.asList(1,2,0,0,0);
        List<Integer> allele11 = Arrays.asList(1,1,0,0,0);
        List<Integer> allele21 = Arrays.asList(2,1,0,0,0);
        List<Integer> allele1 = Arrays.asList(1,0,0,0,0);
        List<Integer> allele2 = Arrays.asList(2,0,0,0,0);
        List<Integer> allele0 = Arrays.asList(0,0,0,0,0);

        //check the ancestral dictionaries
        //todo find a way to extract node numbers in a way that we know their position in the tree
        assertEquals(statesDictionary.get(0).size(),3);
        assertTrue(statesDictionary.get(0).contains(allele12));
        assertTrue(statesDictionary.get(0).contains(allele1));
        assertTrue(statesDictionary.get(0).contains(allele0));

        //2nd leaf
        assertEquals(statesDictionary.get(1).size(),3);
        assertTrue(statesDictionary.get(1).contains(allele11));
        assertTrue(statesDictionary.get(1).contains(allele1));
        assertTrue(statesDictionary.get(1).contains(allele0));

        // node c
        assertEquals(statesDictionary.get(2).size(),3);
        assertTrue(statesDictionary.get(2).contains(allele21));
        assertTrue(statesDictionary.get(2).contains(allele2));
        assertTrue(statesDictionary.get(2).contains(allele0));

        //internal node between a and b
        assertEquals(statesDictionary.get(3).size(),2);
        assertTrue(statesDictionary.get(3).contains(allele1));
        assertTrue(statesDictionary.get(3).contains(allele0));

        //root node a/b/c
        assertEquals(statesDictionary.get(4).size(),1);
        assertTrue(statesDictionary.get(4).contains(allele0));



    }


    @Test
    public void test_likelihood_cherry_homogeneous_identical_edits() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();
        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);


        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);
        RealParameter origin = new RealParameter("6");
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel,"originTime", origin);


        //initialise probabilities
        likelihood.probabilities = new double[tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0,0,0,0,0);
        List<Integer> allele12 = Arrays.asList(1,2,0,0,0);
        List<Integer> allele1 = Arrays.asList(1,0,0,0,0);
        double clock_rate = 0.5;

        double p0000_internal = submodel.getSequenceTransitionProbability(allele0,allele12,5*clock_rate)*submodel.getSequenceTransitionProbability(allele0,allele12,5*clock_rate);
        double p1000_internal = submodel.getSequenceTransitionProbability(allele1,allele12,5*clock_rate)*submodel.getSequenceTransitionProbability(allele1,allele12,5*clock_rate);
        double p1200_internal = submodel.getSequenceTransitionProbability(allele12,allele12,5*clock_rate)*submodel.getSequenceTransitionProbability(allele12,allele12,5*clock_rate);

        //root node
        double proot = p0000_internal * submodel.getSequenceTransitionProbability(allele0,allele0,1*clock_rate) + p1000_internal* submodel.getSequenceTransitionProbability(allele0,allele1,1*clock_rate) + p1200_internal * submodel.getSequenceTransitionProbability(allele0,allele12,1*clock_rate) ;

        //loglikelihood
        double LogPExpected = Math.log(proot);


        double LogPCalc = likelihood.calculateLogP();

        assertEquals(LogPCalc, LogPExpected);


    }

    @Test
    public void test_likelihood_cherry_homogeneous_different_edits() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "(CHILD1:5,CHILD2:5)";
        //Using Origin Time: total height is 6




        Sequence a = new Sequence("CHILD1", "1,2,2,0,0");
        Sequence b = new Sequence("CHILD2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();
        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);
        RealParameter origin = new RealParameter("6");
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel,"originTime", origin);


        //initialise probabilities
        likelihood.probabilities = new double[tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0,0,0,0,0);
        List<Integer> allele12 = Arrays.asList(1,2,0,0,0);
        List<Integer> allele122 = Arrays.asList(1,2,2,0,0);
        List<Integer> allele1 = Arrays.asList(1,0,0,0,0);
        double clock_rate = 0.5;

        double p0000_internal = submodel.getSequenceTransitionProbability(allele0,allele12,5*clock_rate)*submodel.getSequenceTransitionProbability(allele0,allele122,5*clock_rate);
        double p1000_internal = submodel.getSequenceTransitionProbability(allele1,allele12,5*clock_rate)*submodel.getSequenceTransitionProbability(allele1,allele122,5*clock_rate);
        double p1200_internal = submodel.getSequenceTransitionProbability(allele12,allele12,5*clock_rate)*submodel.getSequenceTransitionProbability(allele12,allele122,5*clock_rate);

        //root node
        double proot = p0000_internal * submodel.getSequenceTransitionProbability(allele0,allele0,1*clock_rate) + p1000_internal* submodel.getSequenceTransitionProbability(allele0,allele1,1*clock_rate) + p1200_internal * submodel.getSequenceTransitionProbability(allele0,allele12,1*clock_rate) ;

        //loglikelihood
        double LogPExpected = Math.log(proot);


        double LogPCalc = likelihood.calculateLogP();

        assertEquals(LogPCalc, LogPExpected);


    }

    @Test
    public void test_likelihood_cherry_homogeneous_nothing() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "(CHILD1:5,CHILD2:5)";
        Sequence a = new Sequence("CHILD1", "0,0,0,0,0");
        Sequence b = new Sequence("CHILD2", "0,0,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();
        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);
        RealParameter origin = new RealParameter("6");
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel,"originTime",origin);


        //initialise probabilities
        likelihood.probabilities = new double[tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0,0,0,0,0);
        double clock_rate = 0.5;

        double p0000_internal = submodel.getSequenceTransitionProbability(allele0,allele0,5*clock_rate)*submodel.getSequenceTransitionProbability(allele0,allele0,5*clock_rate);

        //root node
        double proot = p0000_internal * submodel.getSequenceTransitionProbability(allele0,allele0,1*clock_rate) ;

        //loglikelihood
        double LogPExpected = Math.log(proot);


        double LogPCalc = likelihood.calculateLogP();

        assertEquals(LogPCalc, LogPExpected);


    }

    //todo remove/move these tests to testSubModel, as a single branch does not make sense for testing the likelihood
//    @Test
//    public void test_likelihood_cherry_homogeneous_single_branch() {
//
//
//        //Testing the ancestral state reconstruction at internal nodes
//
//        String newick = "(CHILD2:5)";
//        Sequence a = new Sequence("CHILD2", "122000");
//
//        Alignment alignment = new Alignment();
//        alignment.initByName("sequence", a, "dataType", "integer");
//
//        Tree tree1 = new TreeParser();
//        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
//                newick,
//                "adjustTipHeights", false, "offset", 0);
//
//        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();
//
//        //create a sub model with values
//        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();
//        RealParameter freqs = new RealParameter("0.8 0.2");
//        Frequencies frequencies = new Frequencies();
//        frequencies.initByName("frequencies", freqs, "estimate", false);
//        submodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);
//
//        //site model
//        SiteModel siteM = new SiteModel();
//        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);
//
//        //likelihood class
//        RealParameter meanRate = new RealParameter("0.5");
//        StrictClockModel clockModel = new StrictClockModel();
//        RealParameter origin = new RealParameter("5");
//        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel,"originTime",origin);
//
//
//        //initialise probabilities
//        likelihood.probabilities = new double[tree1.getNodeCount()][];
//
//        //Manually calc the likelihood for that tree:
//
//        //internal node partials:
//        List<Integer> allele0 = Arrays.asList(0,0,0,0,0);
//        List<Integer> allele122 = Arrays.asList(1,2,2,0,0);
//        double clock_rate = 0.5;
//
//        double p0000_internal = submodel.getSequenceTransitionProbability(allele0,allele122,5*clock_rate);
//
//        //root node
//        double proot = p0000_internal ;
//
//        //loglikelihood
//        double LogPExpected = Math.log(proot);
//
//
//        double LogPCalc = likelihood.calculateLogP();
//
//        assertEquals(LogPCalc, LogPExpected);
//
//
//    }


//    @Test
//    public void test_likelihood_cherry_homogeneous_single_branch_nothing() {
//
//
//        //Testing the ancestral state reconstruction at internal nodes
//        //single branch
//        String newick = "(CHILD2:5);";
//        Sequence a = new Sequence("CHILD2", "00000");
//
//        Alignment alignment = new Alignment();
//        alignment.initByName("sequence", a, "dataType", "integer");
//
//        Tree tree1 = new TreeParser();
//        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
//                newick,
//                "adjustTipHeights", false, "offset", 0);
//
//        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();
//
//        //create a sub model with values
//        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();
//        RealParameter meanRate = new RealParameter("0.5");
//        StrictClockModel clockModel = new StrictClockModel();
//        clockModel.initByName("clock.rate",meanRate);
//
//
//        RealParameter freqs = new RealParameter("0.8 0.2");
//
//        Frequencies frequencies = new Frequencies();
//        frequencies.initByName("frequencies", freqs, "estimate", false);
//        submodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);
//        //site model
//        SiteModel siteM = new SiteModel();
//        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);
//
//        //likelihood class
//        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel);
//
//
//        //initialise probabilities
//        likelihood.probabilities = new double[tree1.getNodeCount()][];
//
//        //Manually calc the likelihood for that tree:
//
//        //internal node partials:
//        List<Integer> allele0 = Arrays.asList(0,0,0,0,0);
//        double clock_rate = 0.5;
//        double p0000_internal = submodel.getSequenceTransitionProbability(allele0,allele0,5*clock_rate);
//
//        //root node
//        double proot = p0000_internal ;
//
//        //loglikelihood
//        double LogPExpected = Math.log(proot);
//
//
//        double LogPCalc = likelihood.calculateLogP();
//
//        assertEquals(LogPCalc, LogPExpected);
//
//
//    }

    @Test
    public void test_likelihood_3_leaves_homogeneous() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:1,CHILD2:2.0)";

        Sequence a = new Sequence("CHILD1", "1,2,0,0,0");
        Sequence b = new Sequence("CHILD3", "1,2,0,0,0");
        Sequence c = new Sequence("CHILD2", "2,1,0,0,0");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence", b,"sequence", c,"dataType", "integer");


        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();
        RealParameter freqs = new RealParameter("0.8 0.2");

        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);
        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);
        RealParameter origin = new RealParameter("4");
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel,"originTime",origin);




        //initialise probabilities
        likelihood.probabilities = new double[tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0,0,0,0,0);
        List<Integer> allele12 = Arrays.asList(1,2,0,0,0);
        List<Integer> allele21 = Arrays.asList(2,1,0,0,0);
        List<Integer> allele1 = Arrays.asList(1,0,0,0,0);
        double clock_rate = 0.5;

        double p0000_internal1 = submodel.getSequenceTransitionProbability(allele0,allele12,1*clock_rate)*submodel.getSequenceTransitionProbability(allele0,allele12,1*clock_rate);
        double p1000_internal1 = submodel.getSequenceTransitionProbability(allele1,allele12,1*clock_rate)*submodel.getSequenceTransitionProbability(allele1,allele12,1*clock_rate);

        double p1200_internal1 = submodel.getSequenceTransitionProbability(allele12,allele12,1*clock_rate)*submodel.getSequenceTransitionProbability(allele12,allele12,1*clock_rate);

        double p0000_internal2 = (p0000_internal1 * submodel.getSequenceTransitionProbability(allele0,allele0,1*clock_rate) + p1000_internal1 * submodel.getSequenceTransitionProbability(allele0,allele1,1*clock_rate) + p1200_internal1 * submodel.getSequenceTransitionProbability(allele0,allele12,1*clock_rate)) * ( submodel.getSequenceTransitionProbability(allele0,allele21,2*clock_rate));

        //root node
        double proot = p0000_internal2 * submodel.getSequenceTransitionProbability(allele0,allele0,2*clock_rate) ;

        //loglikelihood
        double LogPExpected = Math.log(proot);

        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected,LogPCalc);


    }

    @Test
    public void test_likelihood_debugging_NAN() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((137:8.952403402204816E-4,78:8.952403402204816E-4)325:0.0010347828136821708,(297:3.6835646305956847E-4,149:3.6835646305956847E-4)326:0.0015616666908430839)327:0.029335608497724822";

        Sequence a = new Sequence("137", "2,1,2,0,0");
        Sequence b = new Sequence("78", "1,1,0,0,0");
        Sequence d = new Sequence("297", "2,1,1,1,2");
        Sequence e = new Sequence("149", "2,2,1,0,0");


        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence", b,"sequence", d,"sequence", e,"dataType", "integer");


        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();
        RealParameter freqs = new RealParameter("0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1");

        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);
        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class
        RealParameter meanRate = new RealParameter("0.1");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);
        RealParameter origin = new RealParameter("25");
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel,"originTime",origin);




        //initialise probabilities
        likelihood.probabilities = new double[tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0,0,0,0,0);
        List<Integer> allele12 = Arrays.asList(1,2,0,0,0);
        List<Integer> allele21 = Arrays.asList(2,1,0,0,0);
        List<Integer> allele1 = Arrays.asList(1,0,0,0,0);
        double clock_rate = 0.5;

        double p0000_internal1 = submodel.getSequenceTransitionProbability(allele0,allele12,1*clock_rate)*submodel.getSequenceTransitionProbability(allele0,allele12,1*clock_rate);
        double p1000_internal1 = submodel.getSequenceTransitionProbability(allele1,allele12,1*clock_rate)*submodel.getSequenceTransitionProbability(allele1,allele12,1*clock_rate);

        double p1200_internal1 = submodel.getSequenceTransitionProbability(allele12,allele12,1*clock_rate)*submodel.getSequenceTransitionProbability(allele12,allele12,1*clock_rate);

        double p0000_internal2 = (p0000_internal1 * submodel.getSequenceTransitionProbability(allele0,allele0,1*clock_rate) + p1000_internal1 * submodel.getSequenceTransitionProbability(allele0,allele1,1*clock_rate) + p1200_internal1 * submodel.getSequenceTransitionProbability(allele0,allele12,1*clock_rate)) * ( submodel.getSequenceTransitionProbability(allele0,allele21,2*clock_rate));

        //root node
        double proot = p0000_internal2 * submodel.getSequenceTransitionProbability(allele0,allele0,2*clock_rate) ;

        //loglikelihood
        double LogPExpected = Math.log(proot);

        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected,LogPCalc);


    }

    @Test
    public void test_likelihood_3_leaves_all_different() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((CHILD1:1,CHILD3:1)INTERNAL:1,CHILD2:2.0)";

        Sequence a = new Sequence("CHILD1", "1,1,0,0,0");
        Sequence b = new Sequence("CHILD3", "1,2,0,0,0");
        Sequence c = new Sequence("CHILD2", "2,1,0,0,0");
        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence", b,"sequence", c,"dataType", "integer");


        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();
        RealParameter freqs = new RealParameter("0.8 0.2");


        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class
        
        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);
        RealParameter origin = new RealParameter("4");
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel,"originTime",origin);



        //initialise probabilities
        likelihood.probabilities = new double[tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        List<Integer> allele0 = Arrays.asList(0,0,0,0,0);
        List<Integer> allele12 = Arrays.asList(1,2,0,0,0);
        List<Integer> allele11 = Arrays.asList(1,1,0,0,0);
        List<Integer> allele21 = Arrays.asList(2,1,0,0,0);
        List<Integer> allele1 = Arrays.asList(1,0,0,0,0);
        double clock_rate = 0.5;

        double p0000_internal1 = submodel.getSequenceTransitionProbability(allele0,allele12,1*clock_rate)*submodel.getSequenceTransitionProbability(allele0,allele11,1*clock_rate);
        double p1000_internal1 = submodel.getSequenceTransitionProbability(allele1,allele12,1*clock_rate)*submodel.getSequenceTransitionProbability(allele1,allele11,1*clock_rate);

        double p0000_internal2 = (p0000_internal1 * submodel.getSequenceTransitionProbability(allele0,allele0,1*clock_rate) + p1000_internal1 * submodel.getSequenceTransitionProbability(allele0,allele1,1*clock_rate )) * ( submodel.getSequenceTransitionProbability(allele0,allele21,2*clock_rate));


        //root node
        double proot = p0000_internal2 * submodel.getSequenceTransitionProbability(allele0,allele0,2*clock_rate) ;

        //loglikelihood
        double LogPExpected = Math.log(proot);

        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPExpected,LogPCalc);


    }

    @Test
    public void test_likelihood_debug() {


        //Testing the ancestral state reconstruction at internal nodes
        //tree with 3 tips
        String newick = "((1:0.08286956714231264,0:0.08286956714231264)3:0.23559810414018906,2:0.3184676712825017)4:0.0";

        Sequence a = new Sequence("0", "2,0,0,0,0");
        Sequence b = new Sequence("1", "2,2,1,0,0");
        Sequence c = new Sequence("2", "2,7,1,0,0");



        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a,"sequence", b,"sequence", c,"dataType", "integer");


        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModelHomogeneous submodel = new TypewriterSubstitutionModelHomogeneous();
        RealParameter freqs = new RealParameter("0.09839173 0.09471644 0.08676738 0.07392497 0.05575022 0.04823091 0.03505629 0.02870759 0.02532363 0.02211764 0.01668117 0.01032324 0.004736481 0.004293548 0.004135358 0.003263993 0.003008252 0.00199188 0.0007105381 0.3818687548");


        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class

        RealParameter meanRate = new RealParameter("0.5");
        StrictClockModel clockModel = new StrictClockModel();
        clockModel.initByName("clock.rate",meanRate);
        RealParameter origin = new RealParameter("4");
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel,"originTime",origin);



        //initialise probabilities
        likelihood.probabilities = new double[tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
//        List<Integer> allele0 = Arrays.asList(0,0,0,0,0);
//        List<Integer> allele12 = Arrays.asList(1,2,0,0,0);
//        List<Integer> allele11 = Arrays.asList(1,1,0,0,0);
//        List<Integer> allele21 = Arrays.asList(2,1,0,0,0);
//        List<Integer> allele1 = Arrays.asList(1,0,0,0,0);
//        double clock_rate = 0.5;
//
//        double p0000_internal1 = submodel.getSequenceTransitionProbability(allele0,allele12,1*clock_rate)*submodel.getSequenceTransitionProbability(allele0,allele11,1*clock_rate);
//        double p1000_internal1 = submodel.getSequenceTransitionProbability(allele1,allele12,1*clock_rate)*submodel.getSequenceTransitionProbability(allele1,allele11,1*clock_rate);
//
//        double p0000_internal2 = (p0000_internal1 * submodel.getSequenceTransitionProbability(allele0,allele0,1*clock_rate) + p1000_internal1 * submodel.getSequenceTransitionProbability(allele0,allele1,1*clock_rate )) * ( submodel.getSequenceTransitionProbability(allele0,allele21,2*clock_rate));
//
//
//        //root node
//        double proot = p0000_internal2 * submodel.getSequenceTransitionProbability(allele0,allele0,2*clock_rate) ;
//
//        //loglikelihood
//        double LogPExpected = Math.log(proot);

        double LogPCalc = likelihood.calculateLogP();
//        assertEquals(LogPExpected,LogPCalc);


    }

    ////OLD SUBSTIRUTION MODEL IMPLEMENTATION
//    @Test
//    public void test_likelihood_cherry() {
//
//
//        // Testing the ancestral state reconstruction at internal nodes
//        //tree with 2 tips
//        String newick = "(CHILD1:5,CHILD2:5)";
//        Sequence a = new Sequence("CHILD1", "0102030000");
//        Sequence b = new Sequence("CHILD2", "0102000000");
//
//        Alignment alignment = new Alignment();
//        alignment.initByName("sequence", a, "dataType", "integer");
//        alignment.initByName("sequence", b, "dataType", "integer");
//
//        Tree tree1 = new TreeParser();
//        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
//                newick,
//                "adjustTipHeights", false, "offset", 0);
//
//        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();
//
//        //create a sub model with values
//        TypewriterSubstitutionModel submodel = new TypewriterSubstitutionModel();
//        RealParameter insertrates = new RealParameter("0.5 0.5 0.5 0.0");
//        RealParameter freqs = new RealParameter("1.0 0 0 0 0");
//        Frequencies frequencies = new Frequencies();
//        frequencies.initByName("frequencies", freqs, "estimate", false);
//        submodel.initByName("rates", insertrates, "frequencies", frequencies);
//        submodel.calculateIntermediates();
//
//        //site model
//        SiteModel siteM = new SiteModel();
//        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);
//
//        //likelihood class
//        RealParameter meanRate = new RealParameter("0.5");
//        StrictClockModel clockModel = new StrictClockModel();
//        clockModel.initByName("clock.rate",meanRate);
//        RealParameter origin = new RealParameter("6");
//        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM,"branchRateModel",clockModel,"originTime", origin);
//
//
//
//        //initialise probabilities
//        likelihood.probabilities = new double[tree1.getNodeCount()][];
//
//        //Manually calc the likelihood for that tree:
//
//        //internal node partials:
//        double p0000_internal = submodel.getTransitionProbability(1,5) * submodel.getTransitionProbability(2,5) * submodel.getTransitionProbability(1,5)*submodel.getTransitionProbability(2,5)*submodel.getTransitionProbability(3,5);
//        double p1000_internal = submodel.getTransitionProbability(2,5) * submodel.getTransitionProbability(2,5) * submodel.getTransitionProbability(3,5);
//        double p1200_internal = submodel.getTransitionProbability(0,5) * submodel.getTransitionProbability(3,5);
//
//        //root node
//        double proot = p0000_internal *  submodel.getTransitionProbability(0,1) + p1000_internal* submodel.getTransitionProbability(1,1) + p1200_internal * submodel.getTransitionProbability(1,1) * submodel.getTransitionProbability(2,1);
//
//        //loglikelihood
//        double LogPExpected = Math.log(proot);
//
//        double LogPCalc = likelihood.calculateLogP();
//        assertEquals(LogPCalc, LogPExpected);
//
//
//    }








}

