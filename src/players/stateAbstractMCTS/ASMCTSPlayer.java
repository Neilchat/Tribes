package players.stateAbstractMCTS;

import core.Types;
import core.actions.Action;
import core.actions.tribeactions.EndTurn;
import core.actors.City;
import core.actors.units.Unit;
import core.game.GameState;
import players.Agent;
import players.stateAbstractMCTS.ASMCTSParams;
import players.stateAbstractMCTS.TreeNode;
import utils.ElapsedCpuTimer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class ASMCTSPlayer extends Agent {

    private Random m_rnd;
    private ASMCTSParams params;

    public ASMCTSPlayer(long seed)
    {
        super(seed);
        m_rnd = new Random(seed);
        this.params = new ASMCTSParams();
    }

    public ASMCTSPlayer(long seed, ASMCTSParams params) {
        this(seed);
        this.params = params;
    }

    public Action act(GameState gs, ElapsedCpuTimer ect) {
        //Gather all available actions:
        ArrayList<Action> allActions = gs.getAllAvailableActions();

        if(allActions.size() == 1)
            return allActions.get(0); //EndTurn, it's possible.

//        ArrayList<Action> cityActions = gs.getAllCityActions();
//        ArrayList<Action> unitActions = gs.getAllUnitActions();
//        ArrayList<Action> tribeActions = gs.getTribeActions();
//
//        boolean unitFirst = true;
//        if (m_rnd.nextInt()%2 == 0)
//            tribeActions.addAll(unitActions);
//        else {
//            tribeActions.addAll(cityActions);
//            unitFirst = false;
//        }

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
        if (isClose(opponentUnits, cities)!=null){
            abs = "defend";
        }
        else if (isClose(units, opponentCities)!=null){
            abs = "capture";
        } else if (cityCanUpgrade){
            abs = "upgrade";
        }
        else abs = "advance";

        //todo add kill if enemy near
        //todo if water nearby research port

//        gs.getCities(this.playerID)
//        for (City city :gs.getCities(this.playerID)){
//            city.canLevelUp()
//        }
        TreeNode m_root = new TreeNode(params, m_rnd, allActions.size(), allActions, this.playerID, true, abs, isClose(opponentUnits, cities), isClose(units, opponentCities));
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

    private boolean waterNoDock(GameState gs){
        int count = 0;

        for (int i = 0 ; i< gs.getBoard().getSize(); i++) {
            for (int j = 0; j < gs.getBoard().getSize(); j++) {
                if (gs.getBoard().getTerrainAt(i,j) == Types.TERRAIN.SHALLOW_WATER ||gs.getBoard().getTerrainAt(i,j) == Types.TERRAIN.DEEP_WATER) count ++;
            }
        }
        return true;
    }




    @Override
    public Agent copy() {
        return null;
    }

}