package evolution.datatype;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.*;
import beast.evolution.datatype.DataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.SubstitutionModel;
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
public void test_single_integer(){

    DataType typewriter = new TypewriterData();

    assertEquals("1", typewriter.getCharacter(01));
    assertEquals("11", typewriter.getCharacter(11));

}

@Test
public void test_multiple_integer(){

        DataType typewriter = new TypewriterData();

        assertEquals("1,11", typewriter.getCharacter(0111));
        //typewriter.getCharacter(01);

        //Alignment alignment = new Alignment();
        //alignment.initByName("sequence", a,"dataType", typewriter);

    }

@Test
public void test_typewriter_data() {
    Sequence a = new Sequence("sequence1", "02010100");
    Sequence b = new Sequence("sequence2", "01020000");

    Alignment alignment = new Alignment();
    alignment.initByName("sequence", a,"dataType","TypewriterData");
    alignment.initByName("sequence", b,"dataType","TypewriterData");

    //trying to test how the data is formatted based on this alignment
    //the pattern implementation reo
    Log.info.println("alignment taxon count :"+ alignment.getTaxonCount());
    Log.info.println("alignment as string  " + alignment.toString(true));

    RealParameter insertrates = new RealParameter("0.5 0.5 0.5 0.5");
    TypewriterSubstitutionModel typewritermodel = new TypewriterSubstitutionModel();
    RealParameter freqs = new RealParameter("1.0 0 0 0");
    Frequencies frequencies = new Frequencies();
    frequencies.initByName("frequencies", freqs,
            "estimate", false);
    typewritermodel.initByName("rates", insertrates,"frequencies", frequencies);
    //attempt at testing the transition probabilities
    double start_time = 2;
    double end_time = 1;
    double branch_rate = 0.5;
    Double probability = typewritermodel.getTransitionProbability(1,branch_rate,start_time,end_time);
    assertEquals(probability,0.15803013970713942, 1e-5);

    //initial tests for getting the ancestral states

    //array representation of sequence 1:
    Log.info.println(alignment.getCounts());
    List<Integer> sequence_a = alignment.getCounts().get(0);
    List<Integer> sequence_b = alignment.getCounts().get(1);
    List<List<Integer>> ancs_sequence_a = get_possible_ancestors(sequence_a);
    List<List<Integer>> ancs_sequence_b = get_possible_ancestors(sequence_b);

    // For a sequence with n sites, there are n_edited_sites + 1 ancestral sequences (all edited position + fully unedited))
    // our toy sequence has 3 edited sites, so it should have 4 possible ancestor states.
    assertEquals(ancs_sequence_a.size(),4, 1e-5);

    // our toy sequence has 3 edited sites, so it should have 4 possible ancestor states.
    //Log.info.println("ancestral states of seq a",ancs_sequence_a);
    Log.info.println("Sequence a ancestral sequences : " + ancs_sequence_a);
    Log.info.println("Sequence b ancestral sequences : " + ancs_sequence_b);

    //attempt at the intersection between the two sets:
    ancs_sequence_b.retainAll(ancs_sequence_a);
    Log.info.println("anc sequences b inter anc sequences a : " + ancs_sequence_b );




}


public List<List<Integer>> get_possible_ancestors(List<Integer> sequence) {
    // to get all possible ancestors we just remove edits 1 by 1 along the barcode, starting from the last edited position
    List<List<Integer>> ancestors = new ArrayList();
    //adding the sequence itself as an ancestor
    ancestors.add(sequence);
    List<Integer> ancestor = new ArrayList<>(sequence);
    for(int i = sequence.size()-1;i >= 0; --i) {
        if(sequence.get(i) != 0) {
            //append possible ancestor list
            ancestor.set(i,0);
            ancestors.add(new ArrayList<>(ancestor));
        }
    }
    return ancestors;
}


}
