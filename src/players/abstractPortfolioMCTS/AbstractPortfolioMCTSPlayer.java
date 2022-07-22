package players.abstractPortfolioMCTS;

import core.Types;
import core.actions.Action;
import core.actors.Building;
import core.actors.City;
import core.actors.units.Unit;
import core.game.GameState;
import players.Agent;
import players.portfolio.ActionAssignment;
import utils.ElapsedCpuTimer;
import utils.Vector2d;
import utils.stats.AIStats;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;


public class AbstractPortfolioMCTSPlayer extends Agent {

    private final Random m_rnd;
    private AbstractPortfolioMCTSParams params;
    private AbstractPortfolioTreeNode m_root;
    private AIStats aiStats;

    public AbstractPortfolioMCTSPlayer(long seed)
    {
        super(seed);
        m_rnd = new Random(seed);
        this.params = new AbstractPortfolioMCTSParams();
        this.aiStats = new AIStats(this.playerID);
    }

    public AbstractPortfolioMCTSPlayer(long seed, AbstractPortfolioMCTSParams params) {
        this(seed);
        this.params = params;
    }

    public Action act(GameState gs, ElapsedCpuTimer ect) {
        //Gather all available actions:
        ArrayList<Action> allActions = gs.getAllAvailableActions();

        if(allActions.size() == 1)
            return allActions.get(0); //EndTurn, it's possible.

//        ArrayList<Action> rootActions = params.PRIORITIZE_ROOT ? determineActionGroup(gs, m_rnd) : allActions;
//        if(rootActions == null)
//            return new EndTurn();

        String abs;
        ArrayList<Unit> opponentUnits = gs.getUnits((this.playerID+1)%2);
        ArrayList<Unit> units = gs.getUnits(this.playerID);
        ArrayList<City> opponentCities = gs.getCities((this.playerID+1)%2);
        ArrayList<City> cities = gs.getCities(this.playerID);

        boolean cityCanUpgrade = false;
        for (City city : cities){
            if (city.getPopulation() >= city.getPopulation_need() - 1) {
                cityCanUpgrade = true;
                break;
            }
        }
        if (nextToVillage(gs, units)!=null) {
            abs = "village";
        }
        else if (isClose(opponentUnits, cities)!=null){
            abs = "defend";
        }
        else if (isClose(units, opponentCities)!=null){
            abs = "capture";
        } else if (cityCanUpgrade){
            abs = "upgrade";
        } else if (waterNoDock(gs)) abs = "dock";
        else abs = "advance";

        m_root = new AbstractPortfolioTreeNode(params, m_rnd, this.playerID, abs, isClose(opponentUnits, cities), isClose(units, opponentCities), nextToVillage(gs,units));
        m_root.setRootGameState(m_root, gs, allPlayerIDs);
        m_root.mctsSearch(ect);

        ActionAssignment act = m_root.getActions().get(m_root.bestAction());
//        ActionAssignment act = m_root.getActions().get(m_root.mostVisitedAction());

        this.updateBranchingFactor(gs);

        return act.getAction();
    }

    private City isClose(ArrayList<Unit> opponentUnits, ArrayList<City> cities){
        for (City city : cities) {
            int attackNum = 0;
            for (Unit unit : opponentUnits) {
                if (unit.getPosition().dist(city.getPosition())<2)
                    attackNum++;
            }
            if (attackNum>0){
                return city;
            }
        }
        return null;
    }

    private Vector2d nextToVillage(GameState gs, ArrayList<Unit> units){
        for (Unit unit: units){
            for (Vector2d pos : unit.getPosition().neighborhood(1, 0, gs.getBoard().getSize())){
                if (gs.getBoard().getTerrainAt(pos.x, pos.y)== Types.TERRAIN.VILLAGE)
                    return pos;
            }
        }
        return null;
    }

    private boolean waterNoDock(GameState gs){
        double count = 0;

        for (int i = 0 ; i< gs.getBoard().getSize(); i++) {
            for (int j = 0; j < gs.getBoard().getSize(); j++) {
                if (gs.getBoard().getTerrainAt(i,j) == Types.TERRAIN.SHALLOW_WATER ||gs.getBoard().getTerrainAt(i,j) == Types.TERRAIN.DEEP_WATER) count ++;
            }
        }
        boolean hasDock = false;
        for (City city : gs.getCities(this.playerID)){
            for (Building building : city.getBuildings()){
                if (building.type== Types.BUILDING.PORT) hasDock = true;
            }
        }
        return (count/(gs.getBoard().getSize()*gs.getBoard().getSize())> 0.4 && gs.getTribeTechTree(this.playerID).isResearched(Types.TECHNOLOGY.SAILING) && !hasDock);
    }

    @Override
    public Agent copy() {
        return null;
    }

    /**
     * Returns the number of actions available for each of the actors, from the perspective of this agent.
     * By default, it's the same as the game state says - but overriding this function allows for pruning analysis.
     */
    public ArrayList<Integer> actionsPerUnit(GameState gs)
    {
        HashMap<Integer, Integer> actionsPerUnit = new HashMap<>();
        for(ActionAssignment aas : m_root.getActions())
        {
            int actorID = aas.getActor().getActorId();

            int n = actionsPerUnit.containsKey(actorID) ? actionsPerUnit.get(actorID) + 1 : 1;
            actionsPerUnit.put(actorID, n);
        }

        ArrayList<Integer> actionCounts = new ArrayList<>();
        for(int id : actionsPerUnit.keySet())
        {
            actionCounts.add(actionsPerUnit.get(id));
        }

        return actionCounts;
    }

    /**
     * Returns the total number of actions available in a game state, from the perspective of this agent.
     * By default, it's the same as the game state says - but overriding this function allows for pruning analysis.
     */
    public int actionsPerGameState(GameState gs)
    {
        return m_root.getActions().size();
    }

    public AbstractPortfolioMCTSParams getParams() {
        return params;
    }

    private void updateBranchingFactor(GameState gameState) {
        ArrayList<Integer> actionCounts = actionsPerUnit(gameState);
        aiStats.addBranchingFactor(gameState.getTick(), actionCounts);
        aiStats.addActionsPerStep(gameState.getTick(), actionsPerGameState(gameState));
    }


    /**
     * Function called at the end of the game. May be used by agents for final analysis.
     * @param reward - final reward for this agent.
     */
    public void result(GameState gs, double reward)
    {
        aiStats.print();
    }

}