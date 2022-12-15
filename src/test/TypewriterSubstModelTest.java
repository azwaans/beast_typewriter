package test;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.substitutionmodel.Frequencies;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;
import lineageTree.substitutionmodel.TypewriterSubstitutionModelHomogeneous;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;


@Description("Test substitution model")
public class TypewriterSubstModelTest {

    @Test
    public void testTransitionProbabilities(){

        // Arrange
        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        typewritermodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);

        Sequence a = new Sequence("cell1", "2,1,0,0,0");
        Sequence b = new Sequence("cell2", "1,2,2,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        // expectedProbability : draw 1 event on a Poisson process bounded to 3
        // 0.60653
        // event 2 with probability 0.2

        Double expectedProbability = 0.36788 * 0.2 ;

        // Act
        Double actualProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,1);

        // Assert
        assertEquals(expectedProbability, actualProbability,0.00001);

    }


    @Test
    public void testTransitionProbabilitiesHomogeneousno_edit(){

        // Arrange
        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();


        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        typewritermodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);


        Sequence a = new Sequence("cell1", "2,1,0,0,0");
        Sequence b = new Sequence("cell2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);


        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5);

        // expectedProbability : draw 0 event on a Poisson process bounded to 3
        // 0.60653

        Double expectedProbability = 0.60653;

        //
        assertEquals(expectedProbability, calculatedProbability, 0.00001);


    }

    @Test
    public void testTransitionProbabilitiesHomogeneoussingle_edit(){

        // Arrange
        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        typewritermodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);


        Sequence a = new Sequence("cell1", "2,1,0,0,0");
        Sequence b = new Sequence("cell2", "1,2,1,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);


        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5);

        // expectedProbability : draw 1 event on a Poisson process bounded to 3, with event frequency 0.8
        // 0.30327 * 0.8
        // 0.242616
        Double expectedProbability = 0.24261;

        //
        assertEquals(expectedProbability, calculatedProbability, 0.00001);


    }


    @Test
    public void testTransitionProbabilitiesHomogeneousMulti_edits(){

        // Arrange
        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        typewritermodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);


        Sequence a = new Sequence("cell1", "2,1,0,0,0");
        Sequence b = new Sequence("cell2", "1,2,1,2,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);


        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5);

        // expectedProbability : draw 1 event on a Poisson process bounded to 3, with event frequency 0.8
        // 0.07582 * 0.8 * 0.2
        // 0.0121312

        Double expectedProbability = 0.0121312;
        assertEquals(expectedProbability, calculatedProbability, 0.00001);


    }

    @Test
    public void testTransitionProbabilitiesHomogeneousSaturation(){

        // Arrange
        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        typewritermodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);


        Sequence a = new Sequence("cell1", "2,1,0,0,0");
        Sequence b = new Sequence("cell2", "1,2,1,2,2");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);

        List<Integer> sequence_b = alignment.getCounts().get(1);
        Log.info.println(sequence_b);


        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5);

        // expectedProbability : draw 3 events on a Poisson process bounded to 3:
        // with event frequency 0.8, 0.2, 0.2
        // Poisson probability: 1 - sum(P(n)) = 1 - (P(0) + P(1) + P(2)) = 1 - (0.60653 + 0.30327 + 0.07582) = 0.01438
        // 0.01438 * (0.2*0.2*0.8)
        // 0.00046016
        Double expectedProbability = 0.00046016;
        assertEquals(expectedProbability, calculatedProbability, 0.00001);


    }

}