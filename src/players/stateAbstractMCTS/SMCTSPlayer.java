package players.stateAbstractMCTS;

import core.Types;
import core.actions.Action;
import core.actors.Building;
import core.actors.City;
import core.actors.units.Unit;
import core.game.GameState;
import players.Agent;
import utils.ElapsedCpuTimer;
import utils.Vector2d;

import java.util.ArrayList;
import java.util.Random;

public class SMCTSPlayer extends Agent {

    private Random m_rnd;
    private SMCTSParams params;

    public SMCTSPlayer(long seed)
    {
        super(seed);
        m_rnd = new Random(seed);
        this.params = new SMCTSParams();
    }

    public SMCTSPlayer(long seed, SMCTSParams params) {
        this(seed);
        this.params = params;
    }

    public Action act(GameState gs, ElapsedCpuTimer ect) {
        ArrayList<Action> allActions = gs.getAllAvailableActions();

        if(allActions.size() == 1)
            return allActions.get(0); //EndTurn, it's possible.

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

        TreeNode m_root = new TreeNode(params, m_rnd, allActions.size(), allActions, this.playerID, abs, isClose(opponentUnits, cities), isClose(units, opponentCities), nextToVillage(gs,units));
        m_root.setRootGameState(m_root, gs, allPlayerIDs);

        m_root.mctsSearch(ect);

        return allActions.get(m_root.mostVisitedAction());

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
                if (gs.getBoard().getTerrainAt(pos.x, pos.y)==Types.TERRAIN.VILLAGE)
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

}