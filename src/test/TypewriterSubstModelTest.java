package test;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;


@Description("Test substitution model")
public class TypewriterSubstModelTest {

    @Test
    public void testTransitionProbabilities(){

        // Arrange
        TypewriterSubstitutionModel typewritermodel = new TypewriterSubstitutionModel();

        RealParameter insertrates = new RealParameter("0.5 0.5 0.5 0.5");
        RealParameter freqs = new RealParameter("1.0 0 0 0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);

        typewritermodel.initByName("rates", insertrates, "frequencies", frequencies);

        double start_time = 2;
        double end_time = 1;
        double branch_rate = 0.5;

        Double expectedProbability = 0.15803013970713942;

        // Act
        Double actualProbability = typewritermodel.getEditTransitionProbability(1, branch_rate, start_time, end_time);

        // Assert
        assertEquals(expectedProbability, actualProbability);

    }

}