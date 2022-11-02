package lineageTree.distributions;


import java.util.*;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;





@Description("tree likelihood for a Typewriter alignment given a generic SiteModel, " +
        "a beast tree and a branch rate model. This is a brute-force approach with no consideration for ")

public class TypewriterTreeLikelihood extends Distribution {

    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

    final public Input<TreeInterface> treeInput = new Input<>("tree", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);

    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    protected SubstitutionModel substitutionModel;
    protected BranchRateModel.Base branchRateModel;
    protected SiteModel.Base m_siteModel;
    protected double[] m_branchLengths;

    public Hashtable<Integer,List<List<Integer>>> ancestralStates ;
    public double[][] probabilities ;



    @Override
    public void initAndValidate() {

        int nodeCount = treeInput.get().getNodeCount();
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();
        branchRateModel = new StrictClockModel();
        m_branchLengths = new double[nodeCount];
        ancestralStates = new Hashtable<>() ;


    }

    public SubstitutionModel getSubstitutionModel() {return substitutionModel;}

    @Override
    public List<String> getArguments() {return null;}

    @Override
    public List<String> getConditions() {return null;}

    @Override
    public void sample(State state, Random random) {}

    @Override
    public double calculateLogP() {

        //1st step : calculate all ancestral states in a Postorder traversal
        // need a data-structure to save these!
        //ugly way to do that: Hashmaps
        final TreeInterface tree = treeInput.get();
        traverseAncestral(tree.getRoot());

        //2nd step : calculate likelihood with these states
        // size of the partial likelihoods at each node = state.

        return 0;

    }


    public void traverseAncestral(Node node) {

        Log.info.println(ancestralStates.size());
        Log.info.println("node ID: " + node.getID() );
        Log.info.println("node number: " + node.getNr() );
        Log.info.println("node is root: " + node.isRoot());



        if(! (node == null) && !node.isRoot()) {


            if (node.isLeaf()) {
                List<List<Integer>> LeafStates = get_possible_ancestors(dataInput.get().getCounts().get(node.getNr()));
                ancestralStates.put(node.getNr(), LeafStates);

            } else {

                final Node child1 = node.getLeft();
                final Node child2 = node.getRight();

                traverseAncestral(child1);
                traverseAncestral(child2);

                List<List<Integer>> AncSetChild1 = ancestralStates.get(child1.getNr());
                List<List<Integer>> AncSetChild2 = ancestralStates.get(child2.getNr());
                List<List<Integer>> AncSetNode = new ArrayList<>(AncSetChild1);
                AncSetNode.retainAll(AncSetChild2);
                ancestralStates.put(node.getNr(), AncSetNode);

            }
        }

        else {

                List<Integer> uneditedState = new ArrayList<Integer>() {{
                    add(0);
                    add(0);
                    add(0);
                    add(0);
                    add(0);
                }};
                List<List<Integer>> unedited = new ArrayList(uneditedState);
                ancestralStates.put(node.getNr(), unedited);
                //this takes care of the stem != root node
                final Node child1 = node.getChild(0);
                traverseAncestral(child1);

        }


    }

    public void traverseLikelihood(Node node) {

    //needs a function to get transition prob between 2 states (that multiplies probabiliyies for multiple sequential
    // editing events
        Log.info.println("probabilities"+ Arrays.deepToString(probabilities));


        if(! (node == null) && !node.isRoot()) {


            if (node.isLeaf()) {

                double[] LeafProbabilities = init_probabilities_leaf(ancestralStates.get(node.getNr()).size());
                probabilities[node.getNr()] = LeafProbabilities;

            } else {

                final Node child1 = node.getLeft();
                final Node child2 = node.getRight();

                traverseLikelihood(child1);
                traverseLikelihood(child2);

                probabilities[node.getNr()] = new double[ancestralStates.get(node.getNr()).size()];

            }
        }

        else {

            //root node!

            //this takes care of the stem != root node
            final Node child1 = node.getChild(0);
            traverseLikelihood(child1);
            double[] probabs = new double[1];
            probabilities[node.getNr()] = probabs;

        }



    }

   public double[] init_probabilities_leaf(int size) {

        double[] leafproba = new double[size];
        Arrays.fill(leafproba,0);
        leafproba[0] = 1;
        return leafproba;
   }


    public static List<List<Integer>> get_possible_ancestors(List<Integer> sequence) {
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
