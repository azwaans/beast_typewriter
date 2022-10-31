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



    //initial tests for getting the ancestral states
    //Integer List representation of sequence 1:
    Log.info.println(alignment.getCounts());
    List<Integer> sequence_a = alignment.getCounts().get(0);
    List<Integer> sequence_b = alignment.getCounts().get(1);
    List<List<Integer>> ancs_sequence_a = TypewriterTreeLikelihood.get_possible_ancestors(sequence_a);
    List<List<Integer>> ancs_sequence_b = TypewriterTreeLikelihood.get_possible_ancestors(sequence_b);

    // For a sequence with n sites, there are n_edited_sites + 1 ancestral sequences (all edited position + fully unedited))
    // our toy sequence a has 3 edited sites, so it should have 4 possible ancestor states.
    assertEquals(ancs_sequence_a.size(), 4, 1e-5);

    //attempt at the intersection between ancestor sets for a and b:
    ancs_sequence_b.retainAll(ancs_sequence_a);
    Log.info.println("anc sequences b inter anc sequences a : " + ancs_sequence_b);
    //that's correct!
}


    @Test
    public void test_ancestors_tree() {


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
        SubstitutionModel submodel = new TypewriterSubstitutionModel();
        SiteModel siteM = new SiteModel();
        siteM.initByName( "gammaCategoryCount", 0, "substModel", submodel);
        likelihood.initByName("data",alignment,"tree",tree1,"siteModel",siteM);
        likelihood.traverseAncestral(tree1.getRoot());
        Hashtable<Integer,List<List<Integer>>> statesDictionary = likelihood.ancestralStates;
        Log.info.println("size ofthe states dictionary" + statesDictionary.size());
        for (int i=0;i<statesDictionary.size();++i) {
            Log.info.println(statesDictionary.get(i));
        }

        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        sequence_a.removeAll(sequence_b);
        Log.info.println("sequence_a minus b"+ sequence_a);

    }


}

