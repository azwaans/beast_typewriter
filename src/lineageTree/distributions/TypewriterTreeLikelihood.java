package lineageTree.distributions;


import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
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


    public void traverseAncestral(Node node){

        if (!node.isLeaf()) {
            final Node child1 = node.getLeft();
            final Node child2 = node.getRight();

            if(node.isRoot()) {
                List<Integer> uneditedState = new ArrayList<Integer>(){{
                    add(0);
                    add(0);
                    add(0);
                    add(0);
                    add(0);
                }};
                List<List<Integer>> unedited = new ArrayList(uneditedState);
                ancestralStates.put(node.getNr(),unedited);
            } else {

                List<List<Integer>> AncSetChild1 = ancestralStates.get(child1.getNr());
                List<List<Integer>> AncSetChild2 = ancestralStates.get(child2.getNr());
                List<List<Integer>> AncSetNode = new ArrayList<>(AncSetChild1);
                AncSetNode.retainAll(AncSetChild2);
                ancestralStates.put(node.getNr(), AncSetNode);
            }

            traverseAncestral(child1);
            traverseAncestral(child2);

        }

        List<List<Integer>> LeafStates = get_possible_ancestors(dataInput.get().getCounts().get(node.getNr()));
        ancestralStates.put(node.getNr(),LeafStates);

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
