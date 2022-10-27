package evolution.datatype;

import beast.core.Description;
import beast.evolution.datatype.DataType.Base;

@Description("Datatype for Typewriter Data")
public class TypewriterData extends Base {

    public TypewriterData() {
        stateCount = 13;
        //for now, we use a integer representation of a Typewriter sequence.
        //todo change this mapping to actual nucleotide insert
        codeLength = 2; // with length 3 nucleotide inserts input: codeLength = 3
        codeMap = "01020304050607080910111213" ; // with length 3 nucleotide inserts input codeMap = ATGCCTTAGTT etc
        //with both codeLength 2 and 3 the resulting representation should be an integer array.
        mapCodeToStateSet = new int[13][];
        for (int i = 0; i < 13; i++) {
            mapCodeToStateSet[i] = new int[1];
            mapCodeToStateSet[i][0] = i;
        }


    }

    @Override
    public String getTypeDescription() {
        return "TypewriterData";
    }
    
    @Override
    public boolean isAmbiguousCode(int code) {
    	return code < 0;
    }

    @Override
    public String getCharacter(int code) {
    	if (code < 0) {
    		return "?";
    	}
    	return code + "";
    }
    
	@Override
	public int[] getStatesForCode(int code) {
		return new int[]{code};
	}
}

