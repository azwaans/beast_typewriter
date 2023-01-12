package test;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.substitutionmodel.Frequencies;
import lineageTree.substitutionmodel.TypewriterSubstitutionModelHomogeneous;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;


@Description("Test substitution model")
public class TypewriterSubstModelTest {

    @Test
    public void testTransitionProbabilities(){

        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        typewritermodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);

        Sequence a = new Sequence("cell1", "1,2,0,0,0");
        Sequence b = new Sequence("cell2", "1,2,2,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        // expectedProbability : draw 1 event on a Poisson process bounded to 3
        // P(1)*0.2

        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(0.5);
        Double expectedProbability = dist.probability(1) * 0.2 ;

        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5);

        // Assert
        assertEquals(expectedProbability, calculatedProbability,0.00001);

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
        Sequence b = new Sequence("cell2", "2,1,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);


        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5);

        // expectedProbability : draw 0 event on a Poisson process bounded to 3
        // P(0)
        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(0.5);
        Double expectedProbability = dist.probability(0);

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
        Sequence b = new Sequence("cell2", "2,1,1,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);


        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5);

        // expectedProbability : draw 1 event on a Poisson process bounded to 3, with event frequency 0.8
        // P(1) * 0.8

        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(0.5);

        Double expectedProbability = dist.probability(1)*0.8;

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
        Sequence b = new Sequence("cell2", "2,1,1,2,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);


        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5);

        // expectedProbability : draw 2 events on a Poisson process bounded to 3, with event frequency 0.8
        // P(2) * 0.8 * 0.2
        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(0.5);
        Double expectedProbability = dist.probability(2) * 0.8 * 0.2;

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
        Sequence b = new Sequence("cell2", "2,1,1,2,2");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);


        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5);
        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(0.5);

        // expectedProbability : draw 3 events on a Poisson process bounded to 3:
        // with event frequency 0.8, 0.2
        // 1 - sum(P(n)) * (0.2*0.2*0.8)

        Double expectedProbability = (1 - (dist.probability(0) + dist.probability(1 ) + dist.probability(2))) * 0.8 * 0.2 * 0.2 ;

        assertEquals(expectedProbability, calculatedProbability, 0.00001);


    }

    @Test
    public void testTransitionProbabilitiesNonsensical(){

        // Arrange
        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter freqs = new RealParameter("0.8 0.2");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        typewritermodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);


        Sequence a = new Sequence("cell1", "2,0,0,0,0");
        Sequence b = new Sequence("cell2", "1,2,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);



        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b,0.5);
        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(0.5);

        // expectedProbability : draw 1 event on a Poisson process bounded to 4:
        // with event frequency 0.8, 0.2
        // P(1) * 0.8

        Double expectedProbability = dist.probability(1)  * 0.8  ;

        assertEquals(expectedProbability, calculatedProbability, 0.00001);
        //THIS IS THE PROBABILITY OF TRANSITION BY ADDING 1, this cannot be!


    }

    @Test
    public void testHighClockNoEdit(){

        // Arrange
        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter freqs = new RealParameter("1.0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        typewritermodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);

        Sequence a = new Sequence("cell1", "0,0,0,0,0");
        Sequence b = new Sequence("cell2", "0,0,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        double clockRate = 40;
        double deltaT = 10;
        double distance = clockRate * deltaT;

        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b, distance);
        Double expectedProbability = 1.91516959671401E-174;

        assertEquals(expectedProbability, calculatedProbability, 1E-16);

    }

    @Test
    public void test0ClockWithEdit(){

        // Arrange
        TypewriterSubstitutionModelHomogeneous typewritermodel = new TypewriterSubstitutionModelHomogeneous();

        RealParameter freqs = new RealParameter("1.0");
        Frequencies frequencies = new Frequencies();
        frequencies.initByName("frequencies", freqs, "estimate", false);
        typewritermodel.initByName( "editfrequencies", freqs, "frequencies" ,frequencies);

        Sequence a = new Sequence("cell1", "0,0,0,0,0");
        Sequence b = new Sequence("cell2", "2,0,0,0,0");

        Alignment alignment = new Alignment();
        alignment.initByName("sequence", a, "dataType", "integer");
        alignment.initByName("sequence", b, "dataType", "integer");

        //internal representation of the sequences for the package:
        List<Integer> sequence_a = alignment.getCounts().get(0);
        List<Integer> sequence_b = alignment.getCounts().get(1);

        double clockRate = 0;
        double deltaT = 10;
        double distance = clockRate * deltaT;

        Double calculatedProbability = typewritermodel.getSequenceTransitionProbability( sequence_a, sequence_b, distance);
        org.apache.commons.math.distribution.PoissonDistribution dist = new PoissonDistributionImpl(distance, );

        Double expectedProbability = 0.0; // dist.probability(0);

        assertEquals(expectedProbability, calculatedProbability, 0.00001);

    }

}