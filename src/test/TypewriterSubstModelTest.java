package test;

import beast.core.Description;
import beast.core.parameter.RealParameter;
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
        TypewriterSubstitutionModel typewritermodel = new TypewriterSubstitutionModel();

        RealParameter insertrates = new RealParameter("0.5 0.5 0.5");
        RealParameter freqs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);

        typewritermodel.initByName("rates", insertrates, "frequencies", frequencies);

        double start_time = 1;
        double end_time = 2;
        double branch_rate = 0.5;

        Double expectedProbability = 0.1758778157529951;

        // Act
        Double actualProbability = typewritermodel.getTransitionProbability( 2, branch_rate*(end_time-start_time));

        // Assert
        assertEquals(expectedProbability, actualProbability);

    }

    @Test
    public void testTransitionProbabilities2single_edit(){

        // Arrange
        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter insertrate = new RealParameter("0.5");
        RealParameter freqs = new RealParameter("0.0 0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);

        typewritermodel.initByName("rate", insertrate, "frequencies", frequencies);

        double start_time = 1;
        double end_time = 2;



        // Act
        Sequence a = new Sequence("cell1", "0201000000");
        Sequence b = new Sequence("cell2", "0102010000");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "TypewriterData");
        alignment.initByName("sequence", b, "dataType", "TypewriterData");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);


        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,(end_time-start_time));

        // expectedProbability : draw 1 event on a Poisson process bounded to 3:, with event frequency 0.8
        // 0.30327 * 0.8
        // 0.242616
        Double expectedProbability = 0.24261;

        //
        assertEquals(expectedProbability, calculatedProbability, 0.00001);


    }


    @Test
    public void testTransitionProbabilities2multi_edits(){

        // Arrange
        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter insertrate = new RealParameter("0.5");
        RealParameter freqs = new RealParameter("0.0 0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);

        typewritermodel.initByName("rate", insertrate, "frequencies", frequencies);

        double start_time = 1;
        double end_time = 2;



        // Act
        Sequence a = new Sequence("cell1", "0201000000");
        Sequence b = new Sequence("cell2", "0102010200");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "TypewriterData");
        alignment.initByName("sequence", b, "dataType", "TypewriterData");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);


        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,(end_time-start_time));

        // expectedProbability : draw 1 event on a Poisson process bounded to 3:, with event frequency 0.8
        // 0.07582 * 0.8 * 0.2
        // 0.0121312
        Double expectedProbability = 0.0121312;

        //
        assertEquals(expectedProbability, calculatedProbability, 0.00001);


    }

}