package test;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.*;
import beast.evolution.datatype.DataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.SubstitutionModel;
import evolution.datatype.TypewriterData;
import lineageTree.distributions.TypewriterTreeLikelihood;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;
import org.junit.Test;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.*;


@Description("Datatype test for Typewriter Data")
public class TypewriterDataTest {

@Test
public void test_single_integer(){

    DataType typewriter = new TypewriterData();

    assertEquals("1", typewriter.getCharacter(01));
    assertEquals("11", typewriter.getCharacter(11));

}

@Test
public void test_typewriter_data() {
    Sequence a = new Sequence("sequence1", "02010100");
    Sequence b = new Sequence("sequence2", "01020000");

    Alignment alignment = new Alignment();
    alignment.initByName("sequence", a, "dataType", "TypewriterData");
    alignment.initByName("sequence", b, "dataType", "TypewriterData");

    //internal representation of the sequences for the package:
    List<Integer> sequence_a = alignment.getCounts().get(0);
    List<Integer> sequence_b = alignment.getCounts().get(1);

    //ancestral sequences
    List<List<Integer>> ancs_sequence_a = TypewriterTreeLikelihood.get_possible_ancestors(sequence_a);
    List<List<Integer>> ancs_sequence_b = TypewriterTreeLikelihood.get_possible_ancestors(sequence_b);

    //manually create ancestral states
    List<Integer> allele21 = Arrays.asList(2, 1, 0, 0);
    List<Integer> allele12 = Arrays.asList(1, 2, 0, 0);
    List<Integer> allele1 = Arrays.asList(1, 0, 0, 0);
    List<Integer> allele2 = Arrays.asList(2, 0, 0, 0);
    List<Integer> allele0 = Arrays.asList(0, 0, 0, 0);

    // For a sequence with n sites, there are n_edited_sites + 1 ancestral sequences (all edited position + fully unedited))
    assertEquals(ancs_sequence_a.size(), 4, 1e-5);

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
    public void test_ancestors() {
        // Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "0102030000");
        Sequence b = new Sequence("CHILD2", "0102000000");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "TypewriterData");
        alignment.initByName("sequence", b, "dataType", "TypewriterData");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModel submodel = new TypewriterSubstitutionModel();
        RealParameter insertrates = new RealParameter("0.5 0.5 0.5 0.0");
        RealParameter freqs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("rates", insertrates, "frequencies", frequencies);
        submodel.calculateIntermediates();

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);

        //test ancestral states sets calculations
        likelihood.traverseAncestral(tree1.getRoot());
        Hashtable<Integer,List<List<Integer>>> statesDictionary = likelihood.ancestralStates;


        //first calculate states dictionary
        //manually create states:
        List<Integer> allele123 = Arrays.asList(1, 2, 3, 0, 0);
        List<Integer> allele12 = Arrays.asList(1, 2, 0, 0, 0);
        List<Integer> allele1 = Arrays.asList(1, 0, 0, 0, 0);
        List<Integer> allele0 = Arrays.asList(0, 0, 0, 0, 0);

        assertEquals(4, statesDictionary.size());

        //1st leaf
        assertTrue(statesDictionary.get(0).contains(allele123));
        assertTrue(statesDictionary.get(0).contains(allele12));
        assertTrue(statesDictionary.get(0).contains(allele1));
        assertTrue(statesDictionary.get(0).contains(allele0));

        //2nd leaf
        assertTrue(statesDictionary.get(1).contains(allele12));
        assertTrue(statesDictionary.get(1).contains(allele1));
        assertTrue(statesDictionary.get(1).contains(allele0));

        //internal node
        assertTrue(statesDictionary.get(2).contains(allele12));
        assertTrue(statesDictionary.get(2).contains(allele1));
        assertTrue(statesDictionary.get(2).contains(allele0));

        //root
        assertEquals(statesDictionary.get(3),allele0);


    }


    @Test
    public void test_likelihood_cherry() {


        // Testing the ancestral state reconstruction at internal nodes
        //tree with 2 tips
        String newick = "((CHILD1:5,CHILD2:5)INTERNAL:1):0.0";
        Sequence a = new Sequence("CHILD1", "0102030000");
        Sequence b = new Sequence("CHILD2", "0102000000");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "TypewriterData");
        alignment.initByName("sequence", b, "dataType", "TypewriterData");

        Tree tree1 = new TreeParser();
        tree1.initByName("IsLabelledNewick", true, "taxa", alignment, "newick",
                newick,
                "adjustTipHeights", false, "offset", 0);

        TypewriterTreeLikelihood likelihood = new TypewriterTreeLikelihood();

        //create a sub model with values
        TypewriterSubstitutionModel submodel = new TypewriterSubstitutionModel();
        RealParameter insertrates = new RealParameter("0.5 0.5 0.5 0.0");
        RealParameter freqs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        submodel.initByName("rates", insertrates, "frequencies", frequencies);
        submodel.calculateIntermediates();

        //site model
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);

        //likelihood class
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);

        //initialise probabilities
        likelihood.probabilities = new double[tree1.getNodeCount()][];

        //Manually calc the likelihood for that tree:

        //internal node partials:
        double p0000_internal = submodel.getTransitionProbability(1,5)* submodel.getTransitionProbability(2,5) * submodel.getTransitionProbability(1,5)*submodel.getTransitionProbability(2,5)*submodel.getTransitionProbability(3,5);
        double p1000_internal = submodel.getTransitionProbability(2,5)* submodel.getTransitionProbability(2,5) * submodel.getTransitionProbability(3,5);
        double p1200_internal = submodel.getTransitionProbability(0,5)*submodel.getTransitionProbability(3,5);

        //root node
        double proot = p0000_internal *  submodel.getTransitionProbability(0,1) + p1000_internal* submodel.getTransitionProbability(1,1) + p1200_internal * submodel.getTransitionProbability(1,1) * submodel.getTransitionProbability(2,1);

        //loglikelihood
        double LogPExpected = Math.log(proot);

        double LogPCalc = likelihood.calculateLogP();
        assertEquals(LogPCalc, LogPExpected);


    }


}

