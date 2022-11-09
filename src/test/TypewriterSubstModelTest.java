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

}