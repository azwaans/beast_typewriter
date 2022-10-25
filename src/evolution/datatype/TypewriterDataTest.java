package evolution.datatype;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.*;
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
public void test_typewriter_data() {
    Sequence a = new Sequence("sequence1", "010201");
    Sequence b = new Sequence("sequence2", "010200");

    Alignment alignment = new Alignment();
    alignment.initByName("sequence", a,"dataType","TypewriterData");
    alignment.initByName("sequence", b,"dataType","TypewriterData");

    RealParameter insertrates = new RealParameter("0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5");
    TypewriterSubstitutionModel typewritermodel = new TypewriterSubstitutionModel();
    RealParameter freqs = new RealParameter("1.0 0 0");
    Frequencies frequencies = new Frequencies();
    frequencies.initByName("frequencies", freqs,
            "estimate", false);
    typewritermodel.initByName("rates", insertrates,"frequencies",frequencies);



}

}
