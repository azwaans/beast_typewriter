package lineageTree.simulation;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.datatype.TypewriterData;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.BEASTClassLoader;
import beast.util.PackageManager;
import beast.util.Randomizer;
import feast.nexus.CharactersBlock;
import feast.nexus.NexusBuilder;
import feast.nexus.TaxaBlock;
import lineageTree.substitutionmodel.TypewriterSubstitutionModel;
import lineageTree.substitutionmodel.TypewriterSubstitutionModelHomogeneous;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


@Description("A more flexible alignment simulator adapted from Tim Vaughan's feast implementation")
public class SimulatedTypeWriterAlignment extends Alignment{

        public Input<Tree> treeInput = new Input<>(
                "tree",
                "Tree down which to simulate sequence evolution.",
                Input.Validate.REQUIRED);

        public Input<SiteModel> siteModelInput = new Input<>(
                "siteModel",
                "Site model to use in simulation.",
                Input.Validate.REQUIRED);

        //TODO rename number of targets
        public Input<Integer> sequenceLengthInput = new Input<>(
                "sequenceLength",
                "Length of sequence to simulate.",
                Input.Validate.REQUIRED);

        public Input<Integer> nrOfInsertionsPerSiteInput = new Input<>(
            "nrOfInsertionsPerSite",
            "Number of insertions to add per site",
            Input.Validate.REQUIRED);

    public Input<RealParameter> originInput = new Input<>(
            "origin", "Start of the process, usually the experiment",
            Input.Validate.OPTIONAL);


        public Input<String> outputFileNameInput = new Input<>(
                "outputFileName",
                "Name of file (if any) simulated alignment should be saved to.");

        private Tree tree;
        private SiteModel siteModel;
        private int seqLength;
        private int nrOfInsertionsPerSite;
        private DataType dataType;
        private double originHeight;

        private String ancestralSeqStr;

        public SimulatedTypeWriterAlignment() {
            sequenceInput.setRule(Input.Validate.OPTIONAL);
        }

        @Override
        public void initAndValidate() {

            tree = treeInput.get();
            siteModel = siteModelInput.get();
            //TODO rename for as above
            seqLength = 1; //TODO rewrite for arbitrary #sites! sequenceLengthInput.get();
            // TODO rename per target
            nrOfInsertionsPerSite = nrOfInsertionsPerSiteInput.get();
            sequences.clear();

            if(originInput.get() != null){
                originHeight = originInput.get().getValue();
            }

            grabDataType();

            simulate();

            super.initAndValidate();

            // Write simulated alignment to disk if required
            if (outputFileNameInput.get() != null) {
                try (PrintStream pstream = new PrintStream(outputFileNameInput.get())) {
                    NexusBuilder nb = new NexusBuilder();
                    nb.append(new TaxaBlock(new TaxonSet(this)));
                    nb.append(new CharactersBlock(this));
                    nb.write(pstream);
                } catch (FileNotFoundException ex) {
                    throw new RuntimeException("Error writing to file "
                            + outputFileNameInput.get() + ".");
                }
            }
        }

        /**
         * Perform actual sequence simulation.
         */
        private void simulate() {
            int nTaxa = tree.getLeafNodeCount();

            TypewriterSubstitutionModelHomogeneous substModel = (TypewriterSubstitutionModelHomogeneous) siteModel.getSubstitutionModel();

            //TODO rename getInsertProbabilities
            double[] transitionProbs = substModel.getInsertionProbs();

            int[][] alignment = new int[nTaxa][nrOfInsertionsPerSite];

            Node root = tree.getRoot();

            //TODO rename root sequence
            int[] parentSequence = new int[nrOfInsertionsPerSite];

            if (originHeight != 0){
                // then parent sequence is sequence at origin and we evolve sequence first down to the root
                double deltaT = originHeight - root.getHeight();
                double clockRate = siteModel.getRateForCategory(0, root);

                int possibleEdits = 5;
                //TODO rename to inserts, i.e nrOfNewInserts
                long nEdits = Randomizer.nextPoisson(deltaT * clockRate);

                int insertionI = 0;
                while (possibleEdits > 0 && nEdits > 0){

                    int newInsertion = Randomizer.randomChoicePDF(transitionProbs) + 1;
                    parentSequence[insertionI] = newInsertion;

                    insertionI++;
                    possibleEdits--;
                    nEdits--;
                }
            }

            //ancestralSeqStr = dataType.encodingToString(parentSequence);

            traverse(root, parentSequence,
                    transitionProbs,
                    alignment);

            for (int leafIdx=0; leafIdx<nTaxa; leafIdx++) {
                String seqString = dataType.encodingToString(alignment[leafIdx]);

                String taxonName;
                if (tree.getNode(leafIdx).getID() != null)
                    taxonName = tree.getNode(leafIdx).getID();
                else
                    taxonName = "t" + leafIdx;

                sequenceInput.setValue(new Sequence(taxonName, seqString), this);
            }
        }

        /**
         * Traverse a tree, simulating a sequence alignment down it.
         *
         * @param node Node of the tree
         * @param parentSequence Sequence at the parent node in the tree
         * @param transitionProbs transition probabilities
         * @param regionAlignment alignment for particular region
         */
        private void traverse(Node node,
                              int[] parentSequence,
                              double[] transitionProbs,
                              int[][] regionAlignment) {


            // ignore categories so far

            for (Node child : node.getChildren()) {

                double deltaT = node.getHeight() - child.getHeight();
                double clockRate = siteModel.getRateForCategory(0, child);

                // Draw characters on child sequence
                int[] childSequence = parentSequence.clone();


                // find site where insertion could happen, i.e. the first '0', otw no more simulation necessary
                int insertionI = 0;
                while ((insertionI < nrOfInsertionsPerSite) && (childSequence[insertionI] !=0)){
                    insertionI++;
                }

                if (insertionI == nrOfInsertionsPerSite){
                    ;
                }else{
                // if parent sequence is not yet fully edited, sample new edit
                    // sample number of new insertions
                    // TODO possibleInserts
                    int possibleEdits = nrOfInsertionsPerSite - insertionI;

                    // TODO numberOfInserts
                    long nEdits = Randomizer.nextPoisson(deltaT * clockRate);

                    while (possibleEdits > 0 && nEdits > 0){

                        int newInsertion = Randomizer.randomChoicePDF(transitionProbs) + 1;
                        childSequence[insertionI] = newInsertion;

                        insertionI++;
                        possibleEdits--;
                        nEdits--;
                    }
                }

                if (child.isLeaf()) {
                    System.arraycopy(childSequence, 0,
                            regionAlignment[child.getNr()], 0, childSequence.length);
                } else {
                    traverse(child, childSequence,
                            transitionProbs,
                            regionAlignment);
                }
            }
        }

        /**
         * HORRIBLE function to identify data type from given description.
         */
        private void grabDataType() {
            if (userDataTypeInput.get() != null) {
                dataType = userDataTypeInput.get();
            } else {

                List<String> dataTypeDescList = new ArrayList<>();
                List<String> classNames = PackageManager.find(beast.evolution.datatype.DataType.class, "beast.evolution.datatype");
                for (String className : classNames) {
                    try {
                        DataType thisDataType = (DataType) BEASTClassLoader.forName(className).newInstance();
                        if (dataTypeInput.get().equals(thisDataType.getTypeDescription())) {
                            dataType = thisDataType;
                            break;
                        }
                        dataTypeDescList.add(thisDataType.getTypeDescription());
                    } catch (ClassNotFoundException
                            | InstantiationException
                            | IllegalAccessException e) {
                    }
                }
                if (dataType == null) {
                    throw new IllegalArgumentException("Data type + '"
                            + dataTypeInput.get()
                            + "' cannot be found.  Choose one of "
                            + Arrays.toString(dataTypeDescList.toArray(new String[0])));
                }
            }
        }
}
