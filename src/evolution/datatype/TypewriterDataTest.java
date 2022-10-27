package evolution.datatype;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.*;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.SubstitutionModel;
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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


@Description("Datatype test for Typewriter Data")


public class TypewriterDataTest {

    @Test
    public void test_typewriter_data() {
        Sequence a = new Sequence("sequence1", "02010100");
        Sequence b = new Sequence("sequence2", "01020000");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "TypewriterData");
        alignment.initByName("sequence", b, "dataType", "TypewriterData");


        RealParameter insertrates = new RealParameter("0.5 0.5 0.5 0.5");
        TypewriterSubstitutionModel typewritermodel = new TypewriterSubstitutionModel();
        RealParameter freqs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs,
                "estimate", false);
        typewritermodel.initByName("rates", insertrates, "frequencies", frequencies);
        //attempt at testing the transition probabilities
        double start_time = 2;
        double end_time = 1;
        double branch_rate = 0.5;

        //todo formally test that
        Double probability = typewritermodel.getTransitionProbability(1, branch_rate, start_time, end_time);
        assertEquals(probability, 0.15803013970713942, 1e-5);

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

    }


}

