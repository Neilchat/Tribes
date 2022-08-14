package players.stateAbstractMCTS;

import core.Types;
import core.actions.Action;
import core.actions.tribeactions.EndTurn;
import core.actors.City;
import core.actors.units.Unit;
import core.game.GameState;
import players.heuristics.StateHeuristic;
import utils.ElapsedCpuTimer;
import utils.Pair;
import utils.Vector2d;

import java.util.*;

import static core.Types.ACTION.END_TURN;

public class TreeNode {

    private ASMCTSParams params;

    private TreeNode root;
    private TreeNode parent;
    private TreeNode[] children;
    private double totValue;
    private int nVisits;
    private Random m_rnd;
    private int m_depth;
    private double[] bounds = new double[]{Double.MAX_VALUE, -Double.MAX_VALUE};
    private int fmCallsCount;
    private int playerID;
    private int absNodeID;
    private String abs;
    private City cityUnderAttack;
    private City cityToAttack;
    private Vector2d villagePos;

    private ArrayList<Action> actions;
    private GameState state;

    private GameState rootState;
    private StateHeuristic rootStateHeuristic;

    //From MCTSPlayer
    TreeNode(ASMCTSParams p, Random rnd, int num_actions, ArrayList<Action> actions, int playerID, String abs, City cityUnderAttack, City cityToAttack, Vector2d villagePos) {
        this(p, null, rnd, num_actions, actions, null, playerID, null, null, abs, cityUnderAttack, cityToAttack, villagePos);
    }

    private TreeNode(ASMCTSParams p, TreeNode parent, Random rnd, int num_actions,
                           ArrayList<Action> actions, StateHeuristic sh, int playerID, TreeNode root, GameState state, String abs, City cityUnderAttack, City cityToAttack, Vector2d villagePos) {
        this.params = p;
        this.fmCallsCount = 0;
        this.parent = parent;
        this.m_rnd = rnd;
        this.actions = actions;
        this.root = root;
        children = new TreeNode[num_actions];
        totValue = 0.0;
        this.playerID = playerID;
        this.state = state;
        if(parent != null) {
            m_depth = parent.m_depth + 1;
            this.rootStateHeuristic = sh;
        }
        else {
            m_depth = 0;
        }
        this.abs = abs;
        this.cityUnderAttack = cityUnderAttack;
        this.cityToAttack = cityToAttack;
        this.villagePos = villagePos;

    }

    void setAbsNodeID(int id){
        this.absNodeID = id;
    }

    void setRootGameState(TreeNode root, GameState gs, ArrayList<Integer> allIDs)
    {
        this.state = gs;
        this.root = root;
        this.rootState = gs;
        this.rootStateHeuristic = params.getStateHeuristic(playerID, allIDs);
    }


    void mctsSearch(ElapsedCpuTimer elapsedTimer) {

        double avgTimeTaken;
        double acumTimeTaken = 0;
        long remaining;
        int numIters = 0;

        int remainingLimit = 5;
        boolean stop = false;

        HashMap<Integer, TreeNode> depthToNode = new HashMap<Integer, TreeNode>();
        HashMap<Integer, List<ArrayList<TreeNode>>> depthToNodeGroups = new HashMap<>();
        HashMap<Integer, ArrayList<TreeNode>> absNodeIDToNodes = new HashMap<>();
        HashMap<Integer, Pair<Double, Integer>> absNodeIDToStats = new HashMap<>();
        HashMap<Integer, Integer> absNodeIDToSize = new HashMap<>();

        absNodeIDToSize.put(0, 1);
        absNodeIDToStats.put(0, new Pair<>(0d,0));
        ArrayList<TreeNode> ng = new ArrayList<>();
        ng.add(this);
        absNodeIDToNodes.put(0, ng);


        depthToNode.put(0, this);

        while(!stop){
//            System.out.println("------- " + root.actions.size() + " -------");
            ElapsedCpuTimer elapsedTimerIteration = new ElapsedCpuTimer();
            TreeNode selected = treePolicy(depthToNode, depthToNodeGroups, absNodeIDToStats, absNodeIDToSize, absNodeIDToNodes);
            double delta = selected.rollOut();
            backUp(selected, delta, absNodeIDToStats, absNodeIDToSize, absNodeIDToNodes);
            numIters++;

            //Stopping condition
            if(params.stop_type == params.STOP_TIME) {
                acumTimeTaken += (elapsedTimerIteration.elapsedMillis()) ;
                avgTimeTaken  = acumTimeTaken/numIters;
                remaining = elapsedTimer.remainingTimeMillis();
                stop = remaining <= 2 * avgTimeTaken || remaining <= remainingLimit;
            }else if(params.stop_type == params.STOP_ITERATIONS) {
                stop = numIters >= params.num_iterations;
            }else if(params.stop_type == params.STOP_FMCALLS)
            {
                stop = fmCallsCount > params.num_fmcalls;
            }
        }
    }

    private TreeNode treePolicy(HashMap depthToNode, HashMap<Integer, List<ArrayList<TreeNode>>> depthToNodeGroups, HashMap<Integer, Pair<Double, Integer>> absNodeIDToStats, HashMap<Integer, Integer> absNodeIDToSize, HashMap<Integer, ArrayList<TreeNode>> absNodeIDToNodes) {

        TreeNode cur = this;

        while (!cur.state.isGameOver() /*&& state.getAllAvailableActions().size() > 1 */ && cur.m_depth < params.ROLLOUT_LENGTH)
        {
            if (cur.notFullyExpanded()) {
                return cur.expand(depthToNode, depthToNodeGroups, absNodeIDToStats,  absNodeIDToSize, absNodeIDToNodes);

            } else {
                cur = cur.uct(absNodeIDToStats, absNodeIDToNodes);
            }
        }

        return cur;
    }

    private int tryForceEnd(GameState state, EndTurn endTurn, int depth)
    {
        boolean willForceEnd = (depth > 0 && (depth % params.FORCE_TURN_END) == 0) && endTurn.isFeasible(state);
        if(!willForceEnd)
            return -1; //Not the time, or not available.

        ArrayList<Action> availableActions = state.getAllAvailableActions();
        int actionIdx = 0;
        while(actionIdx < availableActions.size())
        {
            Action act = availableActions.get(actionIdx);
            if(act.getActionType() == END_TURN)
            {
                //Here's the end turn, return it's index.
                return actionIdx;
            }else actionIdx++;
        }

        //This should not happen, but EndTurn is not available here.
        return -1;
    }

    private TreeNode expand(HashMap depthToNode, HashMap<Integer, List<ArrayList<TreeNode>>> depthToNodeGroups, HashMap<Integer, Pair<Double, Integer>> absNodeIDToStats, HashMap<Integer, Integer> absNodeIDToSize, HashMap<Integer, ArrayList<TreeNode>> absNodeIDToNodes) {

        int bestAction = -1;
        if(bestAction == -1)
        {
            //No turn end, expand
            double bestValue = -1;

            for (int i = 0; i < children.length; i++) {
                double x = m_rnd.nextDouble();
                if (x > bestValue && children[i] == null) {
                    bestAction = i;
                    bestValue = x;
                }
            }
        }

        //Roll the state, create a new node and assign it.
        GameState nextState = state.copy();
        ArrayList<Action> availableActions = getActions(this.m_depth, nextState);
        advance(nextState, availableActions.get(bestAction), true);
        ArrayList<Action> nextActions = getActions(this.m_depth+1, nextState);
        TreeNode tn = new TreeNode(params, this, this.m_rnd, nextActions.size(),
        null, rootStateHeuristic, this.playerID, this.m_depth == 0 ? this : this.root, nextState, this.root.abs, this.root.cityUnderAttack, this.root.cityToAttack, this.root.villagePos);

        updateNodeGroupsAndStats(depthToNode, depthToNodeGroups,  absNodeIDToStats, tn, this.m_depth+1, absNodeIDToSize, absNodeIDToNodes);

        children[bestAction] = tn;
        return tn;
    }

    private void updateNodeGroupsAndStats(HashMap depthToNode, HashMap<Integer, List<ArrayList<TreeNode>>> depthToNodeGroups, HashMap<Integer, Pair<Double, Integer>> absNodeIDToStats, TreeNode tn, int depth, HashMap<Integer, Integer> absNodeIDToSize, HashMap<Integer, ArrayList<TreeNode>> absNodeIDToNodes){
        List<ArrayList<TreeNode>> nodeGroups = depthToNodeGroups.get(depth);
        if (nodeGroups==null) {

            //Make first group
            ArrayList<TreeNode> ng = new ArrayList<>();
            ng.add(tn);
            List<ArrayList<TreeNode>> groups = new ArrayList<ArrayList<TreeNode>>();
            groups.add(ng);
            depthToNodeGroups.put(depth, groups);

            //Update stats
            tn.setAbsNodeID(depth*10000);
            Pair<Double, Integer> stats = new Pair<Double, Integer>(tn.totValue, tn.nVisits);
            absNodeIDToStats.put(depth*10000, stats);
            absNodeIDToSize.put(depth*10000, 1);

            absNodeIDToNodes.put(depth*10000, ng);

        } else {
            boolean groupFound = false;
            //iterate over all abstract nodes (node groups) at that depth
            for (int i =0; i<nodeGroups.size();i++) {
                ArrayList<TreeNode> nodeGroup= nodeGroups.get(i);
                //Check if first element of node group and new node have same parents and are similar using the compare function.
                if (tn.parent.equals(nodeGroup.get(0).parent) && compareStates(nodeGroup.get(0).state, tn.state, depth, tn.abs, tn.parent.state, tn.cityUnderAttack)) {

                    //Add new node to exists abstract node
                    tn.setAbsNodeID(depth*10000+i);
                    nodeGroup.add(tn);

                    //Update Stats of abstract node
                    Pair<Double, Integer> lastStats = absNodeIDToStats.get(depth*10000 + i);
                    Double newTot = lastStats.getFirst()*(nodeGroup.size()-1)/nodeGroup.size() + tn.totValue/nodeGroup.size();
                    Integer newVisit = lastStats.getSecond()*(nodeGroup.size()-1)/nodeGroup.size() + tn.nVisits/nodeGroup.size();
                    absNodeIDToStats.remove(depth*10000 + i);
                    absNodeIDToStats.put(depth*10000 + i, new Pair<>(newTot, newVisit));
                    absNodeIDToSize.put(depth*10000, nodeGroup.size());
                    groupFound = true;
                    break;
                }
            }
            //If no abstract node was found matching
            if (!groupFound) {
                //Create and add new singleton abstract node with the new node
                tn.setAbsNodeID(depth*10000+nodeGroups.size());
                ArrayList<TreeNode> ng = new ArrayList<>();
                ng.add(tn);
                absNodeIDToNodes.put(depth*10000+nodeGroups.size(), ng);
                depthToNodeGroups.get(depth).add(ng);

                //Create the stats for the new node
                Pair<Double, Integer> stats = new Pair<Double, Integer>(tn.totValue, tn.nVisits);
                absNodeIDToStats.put(depth*10000+depthToNodeGroups.get(depth).size()-1, stats);
                absNodeIDToSize.put(depth*10000+depthToNodeGroups.get(depth).size()-1, 1);
            }
        }
    }

    //Similarity Function that compares two states given an abstraction heuristic that can be one of
    // 'defend', 'capture', 'village', 'upgrade', 'dock' or 'advance'
    private boolean compareStates(GameState gs1, GameState gs2, int depth, String abs, GameState parentState, City cityUnderAttack){

        if (abs.equalsIgnoreCase("defend")){

            Defence defParent = getDefenceStats(parentState, cityUnderAttack);
            Defence defGs1 = getDefenceStats(gs1, cityUnderAttack);
            Defence defGs2 = getDefenceStats(gs2, cityUnderAttack);

            boolean isDef1 = false;
            for (Unit unit : gs1.getUnits(this.playerID)){
                if (unit.getPosition().equals(cityUnderAttack.getPosition())) {
                    isDef1 = true;
                    break;
                }
            }
            boolean isDef2 = false;
            for (Unit unit : gs2.getUnits(this.playerID)){
                if (unit.getPosition().equals(cityUnderAttack.getPosition())) {
                    isDef2 = true;
                    break;
                }
            }

            if (parentState.getTribeTechTree(this.playerID).isResearched(Types.TECHNOLOGY.SHIELDS))
                return defGs1.def==defGs2.def;
            else
                return isDef1 == isDef2;
//            }

//            if (depth%4==1){
//                return defGs1.dist==defGs2.dist;
//            }
////
////            if (depth%4==2){
////                return defParent.dist<defGs2.dist;
////            }
//
//            if (depth%4==2){
//                return defGs1.number==defGs2.number;
//            }
//
//            else {
//                return defGs1.health==defGs2.health;
//            }

        } else if (abs.equalsIgnoreCase("capture")){

            Defence defParent = getDefenceStats(parentState, cityToAttack);
            Defence defGs1 = getDefenceStats(gs1, cityToAttack);
            Defence defGs2 = getDefenceStats(gs2, cityToAttack);

//            if (gs1.getCities((this.playerID+1)%2).size() != gs2.getCities((this.playerID+1)%2).size()) {
//                System.out.print(abs);
//            }
//            if (depth%4==0){
//                return defGs1.dist==defGs2.dist;
//            }
//
////            if (defParent.dist<defGs1.dist && defParent.number == defGs1.number){
////                return defParent.dist<defGs2.dist;
////            }
//
//            if (depth%4==1){
//                return defGs1.number==defGs2.number;
//            }
//
//            if (depth%4==2){
//                return defGs1.health!=defGs2.health;
//            }
            boolean isCapturing1 = false;
            boolean isCapturing2 = false;

            for (Unit unit : gs1.getUnits(this.playerID)){
                if (unit.getPosition().equals(cityToAttack.getPosition())) {
                    isCapturing1 = true;
                    break;
                }
            }
            for (Unit unit : gs2.getUnits(this.playerID)){
                if (unit.getPosition().equals(cityToAttack.getPosition())) {
                    isCapturing2 = true;
                    break;
                }
            }
            return isCapturing1 == isCapturing2 || defGs1.dist==defGs2.dist;
        } else if (abs.equalsIgnoreCase("upgrade")){

            int totLevel1 = 0;
            for (City city: gs1.getCities(this.playerID)){
                totLevel1 += city.getLevel();
            }

            int totLevel2 = 0;
            for (City city: gs2.getCities(this.playerID)){
                totLevel2 += city.getLevel();
            }
//            if (totLevel1 != totLevel2) {
//                System.out.print(abs);
//            }
            return totLevel1 == totLevel2;
        } else if (abs.equalsIgnoreCase("village")){
            boolean isCapturing1 = false;
            boolean isCapturing2 = false;

            for (Unit unit : gs1.getUnits(this.playerID)){
                if (unit.getPosition().equals(villagePos)) {
                    isCapturing1 = true;
                    break;
                }
            }
            for (Unit unit : gs2.getUnits(this.playerID)){
                if (unit.getPosition().equals(villagePos)) {
                    isCapturing2 = true;
                    break;
                }
            }
            return isCapturing1 == isCapturing2;

        } else if (abs.equalsIgnoreCase("dock")){
            int count1 = 0;
            for (int i = 0 ; i< gs1.getBoard().getSize(); i++) {
                for (int j = 0; j < gs1.getBoard().getSize(); j++) {
                    if (gs1.getBoard().getBuildingAt(i,j) == Types.BUILDING.PORT) count1 ++;
                }
            }
            int count2 = 0;
            for (int i = 0 ; i< gs2.getBoard().getSize(); i++) {
                for (int j = 0; j < gs2.getBoard().getSize(); j++) {
                    if (gs2.getBoard().getBuildingAt(i,j) == Types.BUILDING.PORT) count2 ++;
                }
            }
            int countPar = 0;
            for (int i = 0 ; i< parentState.getBoard().getSize(); i++) {
                for (int j = 0; j < parentState.getBoard().getSize(); j++) {
                    if (parentState.getBoard().getBuildingAt(i,j) == Types.BUILDING.PORT) countPar ++;
                }
            }
            if (countPar<count1 || countPar<count2)
                return count1==count2;
            return true;
        } else{
//            if (gs1.getTick() < 10) {
            int numRes1 = 0;
            for (Boolean res : gs1.getTribeTechTree(this.playerID).getResearched()) {
                if (res) numRes1 += 1;
            }
            int numRes2 = 0;
            for (Boolean res : gs2.getTribeTechTree(this.playerID).getResearched()) {
                if (res) numRes2 += 1;
            }

//            int rootHealthOther1 = 0;
//            for (Unit unit : parentState.getUnits((this.playerID+1)%2)){
//                rootHealthOther1 += unit.getCurrentHP();
//            }
//
//            int rootHealthOther2 = 0;
//            for (Unit unit : parentState.getUnits((this.playerID+1)%2)){
//                rootHealthOther2 += unit.getCurrentHP();
//            }
//
//            if (depth%3==0) return rootHealthOther1==rootHealthOther2;
//            if (depth%3==1) return gs1.getUnits(this.playerID).size()==gs2.getUnits(this.playerID).size();

            double furthestDist1 = 0;
            int totHealth1 = 0;
            for (City city : gs1.getCities(this.playerID)) {
                if (city.isCapital()) {
                    for (Unit unit : gs1.getUnits(this.playerID)) {
                        totHealth1 +=unit.getCurrentHP();
                        if (unit.getPosition().dist(city.getPosition())> furthestDist1)
                            furthestDist1 = unit.getPosition().dist(city.getPosition());
                    }
                    break;
                }
            }
            double furthestDist2 = 0;
            int totHealth2 = 0;
            for (City city : gs2.getCities(this.playerID )) {
                if (city.isCapital()) {
                    for (Unit unit : gs2.getUnits(this.playerID)) {
                        totHealth2 +=unit.getCurrentHP();
                        if (unit.getPosition().dist(city.getPosition())> furthestDist2)
                            furthestDist2 = unit.getPosition().dist(city.getPosition());
                    }
                    break;
                }
            }

            double distToEnemy1 = 0;
            for (City city : gs1.getCities((this.playerID+1 )%2)) {
                if (city.isCapital()) {
                    for (Unit unit : gs1.getUnits(this.playerID)) {
                        distToEnemy1 = +unit.getPosition().dist(city.getPosition());
                    }
                    break;
                }
            }

            double distToEnemy2 = 0;
            for (City city : gs2.getCities((this.playerID+1 )%2)) {
                if (city.isCapital()) {
                    for (Unit unit : gs2.getUnits(this.playerID)) {
                        distToEnemy2 = +unit.getPosition().dist(city.getPosition());
                    }
                    break;
                }
            }

            int totEnemyHealth1 = 0;
            int totEnemyHealth2 = 0;

            for (Unit unit : gs1.getUnits(this.playerID)) {
                totEnemyHealth1 += unit.getCurrentHP();
            }
            for (Unit unit : gs2.getUnits(this.playerID)) {
                totEnemyHealth2 += unit.getCurrentHP();
            }

            if (gs1.getTick()<10) {
                if (gs1.getTick()%2==0)
                    return numRes1 == numRes2;
                else
                    return furthestDist1 == furthestDist2;
            } else {
                return numRes1 == numRes2;
//                if (gs1.getTick() % 3 == 0)
//                    return distToEnemy1 == distToEnemy2;
//                if (gs1.getTick() % 3 == 1)
//                    return totEnemyHealth1 == totEnemyHealth2;
//                else
//                    return totHealth1 == totHealth2;
            }
        }
    }

    private ArrayList<Action> getActions(int depth, GameState gs){
        return gs.getAllAvailableActions();
    }

    //Gets some stats for a gs, city pair
    private Defence getDefenceStats(GameState gs, City city){

        //Is there a defender at the city
        boolean rootDef = false;
        for (Unit unit : gs.getUnits(this.playerID)){
            if (unit.getType() == Types.UNIT.DEFENDER){
                if (unit.getPosition().equals(city.getPosition())) rootDef = true;
            }
        }

        //Distance of all units from the city
        double rootDist = 0;
        for (Unit unit : gs.getUnits(this.playerID)){
            rootDist += unit.getPosition().dist(city.getPosition());
        }

        //Total number of units
        int rootUnits = gs.getUnits(this.playerID).size();

        //Total enemy health
        int rootHealthOther = 0;
        for (Unit unit : gs.getUnits((this.playerID+1)%2)){
            rootHealthOther += unit.getCurrentHP();
        }
        return new Defence(rootDef, rootDist, rootUnits, rootHealthOther);
    }

    private void advance(GameState gs, Action act, boolean computeActions)
    {
        gs.advance(act, computeActions);
        root.fmCallsCount++;
    }


    private TreeNode uct(HashMap<Integer, Pair<Double, Integer>> absNodeIDToStats, HashMap<Integer, ArrayList<TreeNode>> absNodeIDToNodes) {

        TreeNode selected;
        boolean IamMoving = (state.getActiveTribeID() == this.playerID);


        //No end turn, use uct.
        double[] vals = new double[this.children.length];
        double[] absvals = new double[this.children.length];

        for(int i = 0; i < this.children.length; ++i)
        {
            TreeNode child = children[i];

            //Construct uct values of ground truth nodes, breaks ties randomly
            double hvVal = child.totValue;
            double childValue =  hvVal / (child.nVisits + params.epsilon);
            childValue = normalise(childValue, bounds[0], bounds[1]);
            double uctValue = childValue +
                    params.K * Math.sqrt(Math.log(this.nVisits + 1) / (child.nVisits + params.epsilon));
            uctValue = noise(uctValue, params.epsilon, this.m_rnd.nextDouble());     //break ties randomly

            //Construct uct values of abstract nodes which the children belong to, don't break ties.
            double absHvVal = absNodeIDToStats.get(child.absNodeID).getFirst();
            double absChildValue =  absHvVal / (absNodeIDToStats.get(child.absNodeID).getSecond() + params.epsilon);
            absChildValue = normalise(absChildValue, bounds[0], bounds[1]);
            double absUctValue = absChildValue +
                    params.K * Math.sqrt(Math.log(this.nVisits + 1) / (absNodeIDToStats.get(child.absNodeID).getSecond() + params.epsilon));

            vals[i] = uctValue;
            absvals[i] = absUctValue;
        }

//            if (!Arrays.stream(absvals).max().equals(Arrays.stream(absvals).min())){
//                System.out.println(2);
//            }

        double bestValue = IamMoving ? -Double.MAX_VALUE : Double.MAX_VALUE;
        double bestAbsValue = IamMoving ? -Double.MAX_VALUE : Double.MAX_VALUE;

        //Select best abstract node
        for (double absval : absvals) {
            if ((IamMoving && absval >= bestAbsValue) || (!IamMoving && absval <= bestAbsValue)) {
                bestAbsValue = absval;
            }
        }

        int which = -1;

        //Select best ground truth node inside best abstract node
        for(int i = 0; i < vals.length; ++i) {
            if (absvals[i]!=bestAbsValue) continue;
            if ((IamMoving && vals[i] > bestValue) || (!IamMoving && vals[i] < bestValue)) {
                which = i;
                bestValue = vals[i];
            }
        }


        if (which == -1)
        {
            //if(this.children.length == 0)
            System.out.println("Warning! couldn't find the best UCT value " + which + " : " + this.children.length + " " +
                    //throw new RuntimeException("Warning! couldn't find the best UCT value " + which + " : " + this.children.length + " " +
                    + bounds[0] + " " + bounds[1]);
            System.out.print(this.m_depth + ", AmIMoving? " + IamMoving + ";");
            for(int i = 0; i < this.children.length; ++i)
                System.out.printf(" %f2", vals[i]);
            System.out.println("; selected: " + which);

            which = m_rnd.nextInt(children.length);
        }

        selected = children[which];

//            System.out.print(this.m_depth + ", AmIMoving? " + IamMoving + ";");
//            for(int i = 0; i < this.children.length; ++i)
//                System.out.printf(" %f2", vals[i]);
//            System.out.println("; selected: " + which);



        //Roll the state. This is closed loop, we don't advance the state. We can't do open loop here because the
        // number of actions available on a state depend on the state itself, and random events triggered by multiple
        // runs over the same tree node would have different outcomes (i.e Examine ruins).
        //advance(state, actions.get(selected.childIdx), true);

        root.fmCallsCount++;

        return selected;
    }

    private double rollOut()
    {
        if(params.ROLOUTS_ENABLED) {
            GameState rolloutState = state.copy();
            int thisDepth = this.m_depth;
            while (!finishRollout(rolloutState, thisDepth)) {
                EndTurn endTurn = new EndTurn(rolloutState.getActiveTribeID());
                int bestAction = tryForceEnd(rolloutState, endTurn, thisDepth);
                Action next = (bestAction != -1) ? endTurn : rolloutState.getAllAvailableActions().get(m_rnd.nextInt(rolloutState.getAllAvailableActions().size()));
                advance(rolloutState, next, true);
                thisDepth++;
            }
            return normalise(this.rootStateHeuristic.evaluateState(root.rootState, rolloutState), 0, 1);
        }

        return normalise(this.rootStateHeuristic.evaluateState(root.rootState, this.state), 0, 1);
    }

    private boolean finishRollout(GameState rollerState, int depth)
    {
        if (depth >= params.ROLLOUT_LENGTH)      //rollout end condition.
            return true;

        //end of game
        return rollerState.isGameOver();
    }


    private void backUp(TreeNode node, double result, HashMap<Integer, Pair<Double, Integer>> absNodeIDToStats, HashMap<Integer, Integer> absNodeIDToSize, HashMap<Integer, ArrayList<TreeNode>> absNodeIdToNodes)
    {
        TreeNode n = node;
        while(n != null)
        {
            n.nVisits++;
            n.totValue += result;
            if (result < n.bounds[0]) {
                n.bounds[0] = result;
            }
            if (result > n.bounds[1]) {
                n.bounds[1] = result;
            }

            //Backup abstract nodes
            int size = absNodeIdToNodes.get(n.absNodeID).size();
            Pair<Double, Integer> lastStats = absNodeIDToStats.get(n.absNodeID);
            Double newTot = lastStats.getFirst()*(size-1)/size + n.totValue/size;
            Integer newVisit = lastStats.getSecond()*(size-1)/size + n.nVisits/size;
            absNodeIDToStats.remove(n.absNodeID);
            absNodeIDToStats.put(n.absNodeID, new Pair<>(newTot, newVisit));

            n = n.parent;
        }
    }


    int mostVisitedAction() {
        int selected = -1;
        double bestValue = -Double.MAX_VALUE;
        boolean allEqual = true;
        double first = -1;

        for (int i=0; i<children.length; i++) {

            if(children[i] != null)
            {
                if(first == -1)
                    first = children[i].nVisits;
                else if(first != children[i].nVisits)
                {
                    allEqual = false;
                }

                double childValue = children[i].nVisits;
                childValue = noise(childValue, params.epsilon, this.m_rnd.nextDouble());     //break ties randomly
                if (childValue > bestValue) {
                    bestValue = childValue;
                    selected = i;
                }
            }
        }

        if (selected == -1)
        {
            selected = 0;
        }else if(allEqual)
        {
            //If all are equal, we opt to choose for the one with the best Q.
            selected = bestAction();
        }

        return selected;
    }

    private int bestAction()
    {
        int selected = -1;
        double bestValue = -Double.MAX_VALUE;

        for (int i=0; i<children.length; i++) {

            if(children[i] != null) {
                double childValue = children[i].totValue / (children[i].nVisits + params.epsilon);
                childValue = noise(childValue, params.epsilon, this.m_rnd.nextDouble());     //break ties randomly
                if (childValue > bestValue) {
                    bestValue = childValue;
                    selected = i;
                }
            }
        }

        if (selected == -1)
        {
            System.out.println("Unexpected selection!");
            selected = 0;
        }

        return selected;
    }


    private boolean notFullyExpanded() {
        for (TreeNode tn : children) {
            if (tn == null) {
                return true;
            }
        }

        return false;
    }

    private double normalise(double a_value, double a_min, double a_max)
    {
        if(a_min < a_max)
            return (a_value - a_min)/(a_max - a_min);
        else    // if bounds are invalid, then return same value
            return a_value;
    }

    private double noise(double input, double epsilon, double random)
    {
        return (input + epsilon) * (1.0 + epsilon * (random - 0.5));
    }

    public class Defence{

        public final Boolean def;
        public final Double dist;
        public final int number;
        public final int health;

        public Defence(Boolean def, Double dist, int number, int health){

            this.def = def;
            this.dist = dist;
            this.number = number;
            this.health = health;
        }

    }

}
