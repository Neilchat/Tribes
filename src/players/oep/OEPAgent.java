package players.oep;


import core.actions.Action;
import core.actions.cityactions.CityAction;
import core.actions.unitactions.UnitAction;
import core.actors.Actor;
import core.actors.City;
import core.actors.units.Unit;
import core.game.GameState;
import players.Agent;
import players.heuristics.StateHeuristic;
import utils.ElapsedCpuTimer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class OEPAgent extends Agent {

    private Random m_rnd;
    private StateHeuristic heuristic;
    private OEPParams params;


    private Individual bestIndividual;
    private boolean returnAction = false;
    private int fmCallsCount;

    public OEPAgent(long seed, OEPParams params) {
        super(seed);
        m_rnd = new Random(seed);
        this.params = params;
    }

    @Override
    public Action act(GameState gs, ElapsedCpuTimer ect) {
        double avgTimeTaken;
        double acumTimeTaken = 0;
        long remaining;
        int numIters = 0;

        int remainingLimit = 5;
        boolean stop = false;

        if(this.returnAction){
            Action action;
            if(bestIndividual.getActions().size() == 1){
                returnAction = false;
            }
            action = bestIndividual.returnNext();
            return action;
        }


        //create a population of individuals defined in param
        ArrayList<Individual> population = new ArrayList<>();
        for(int i = 0; i < params.POP_SIZE; i++){
            population.add(new Individual(randomActions(gs.copy())));
        }

        this.heuristic = params.getHeuristic(playerID, allPlayerIDs);

        //keep going until time limit gone
        while(!stop){
            numIters ++;
            ElapsedCpuTimer elapsedTimerIteration = new ElapsedCpuTimer();

            //perform uniform crossover
            boolean even = false;
            if(population.size() % 2 == 0){
                even = true;
            }
            Individual person1 = null;
            if(!even){
                person1 = population.get(m_rnd.nextInt(population.size()));
                population.remove(person1);
            }

            population = procreate(gs.copy(), population);

            if(!even){
                population.add(crossover(gs.copy(), person1, population.get(m_rnd.nextInt(population.size()))));
            }


            // rate each individual and sort them
            for(Individual individual : population){
                individual.setValue(eval(gs.copy(), individual));
            }
            Collections.sort(population, Collections.reverseOrder());
            this.bestIndividual = population.get(population.size() - 1);

            //Kill the amount of the population that needs to die
            int amount =  (int)(population.size() * params.KILL_RATE);
            for(int i = 0; i < amount; i++){
                population.remove(population.size() - 1);
            }

            // fill the population with random individuals
            for(int i = population.size() - 1; i < params.POP_SIZE; i++){
                population.add(new Individual(randomActions(gs.copy())));
            }


            if(params.stop_type == params.STOP_TIME) {
                acumTimeTaken += (elapsedTimerIteration.elapsedMillis()) ;
                avgTimeTaken  = acumTimeTaken/numIters;
                remaining = ect.remainingTimeMillis();
                stop = remaining <= 2 * avgTimeTaken || remaining <= remainingLimit;
                if(stop){ returnAction = true; }
            }else if(params.stop_type == params.STOP_ITERATIONS) {
                stop = numIters >= params.num_iterations;
                if(stop){ returnAction = true; }
            }else if(params.stop_type == params.STOP_FMCALLS){
                stop = fmCallsCount > params.num_fmcalls;
                if(stop){ returnAction = true; }
            }

        }

        Action action = null;
        if(this.returnAction){

            if(bestIndividual.getActions().size() == 1){
                returnAction = false;
            }
            action = bestIndividual.returnNext();

        }

        return action;


    }

    @Override
    public Agent copy() {
        return null;
    }

    private ArrayList<Action> randomActions(GameState gs){
        ArrayList<Action> individual = new ArrayList<>();
        while (!gs.isGameOver() && (gs.getActiveTribeID() == getPlayerID())){
            ArrayList<Action> allAvailableActions = this.allGoodActions(gs, m_rnd);
            Action a = allAvailableActions.get(m_rnd.nextInt(allAvailableActions.size()));
            advance(gs, a);
            individual.add(a);

        }
        return individual;
    }

    // must be given a even no. of population
    private ArrayList<Individual> procreate (GameState gs, ArrayList<Individual> population){
        ArrayList<Individual> group1 = new ArrayList<>();
        ArrayList<Individual> group2 = new ArrayList<>();
        boolean g1 = true;
        while(population.size() > 0){
            Individual temp = population.get(m_rnd.nextInt(population.size()));
            if(g1){
                group1.add(temp);
                g1 = false;
            }else{
                group2.add(temp);
                g1 = true;
            }
            population.remove(temp);
        }

        for(int i = 0; i < group1.size(); i++){
            population.add(crossover(gs.copy(), group1.get(i), group2.get(i)));
        }
        return population;
    }

    //method to perform uniform crossover on two individuals
    private Individual crossover(GameState clone, Individual individual1, Individual individual2){
        ArrayList<Action> in1 = individual1.getActions();
        ArrayList<Action> in2 = individual2.getActions();

        ArrayList<Action> child = new ArrayList<>();
        //if both individuals are of the same size

        boolean ind1 = true;
        boolean sameSize = false;
        int smallSize = in1.size();
        if(smallSize > in2.size()){
            ArrayList<Action> temp = in1;
            in1 = in2;
            in2 = temp;
            smallSize = in1.size();
        }
        else if(smallSize == in2.size()){
            sameSize = true;
        }
        smallSize --;
        if(!sameSize && in1.size() > 0){in1.remove(in1.size() - 1);}
        int in1amount =(int)(in1.size() / 2);
        int in2amount = in1.size() - in1amount;
        for(int i = 0; i < in2.size(); i++){
            if(i >= smallSize){
                child.add(in2.get(i));
            }else{
                if(in1amount == 0){
                    child.add(in2.get(i));
                }else if(in2amount == 0){
                    child.add(in1.get(i));
                }else{
                    int temp = m_rnd.nextInt(100);
                    if(temp < 50){
                        child.add(in1.get(i));
                        in1amount--;
                    }else{
                        child.add(in2.get(i));
                        in2amount--;
                    }
                }
            }
        }
        child = repair(clone, child);
        return (new Individual(child));
    }
    //repair an individual if actions can't be performed with a random action
    private ArrayList<Action> repair(GameState gs, ArrayList<Action> child){
        ArrayList<Action> repairedChild = new ArrayList<>();
        boolean mutated = false;
        for(int a = 0 ;a < child.size(); a ++) {
            if(!(gs.getActiveTribeID() == getPlayerID())){return repairedChild;}
            int chance = m_rnd.nextInt((int)(params.MUTATION_RATE * 100));
            if((m_rnd.nextInt(100) < chance) && !mutated){
                Action ac = mutation(gs);
                advance(gs, ac);
                repairedChild.add(ac);
                mutated = true;
            }else{
                GameState copy = gs.copy();
                boolean added = false;
                try {
                    boolean done = checkActionFeasibility(child.get(a), gs.copy());
                    if (!done) {
                        ArrayList<Action> allAvailableActions = this.allGoodActions(gs.copy(), m_rnd);
                        Action ac = allAvailableActions.get(m_rnd.nextInt(allAvailableActions.size()));
                        advance(gs,ac);
                        repairedChild.add(ac);
                        added = true;
                    } else {
                        repairedChild.add(child.get(a));
                        added = true;
                        advance(gs, child.get(a));
                    }
                } catch (Exception e) {
                    if(added){repairedChild.remove(repairedChild.size()-1);}
                    gs = copy;

//                    ArrayList<Action> allAvailableActions = this.allGoodActions(gs, m_rnd);
//                    Action ac = allAvailableActions.get(m_rnd.nextInt(allAvailableActions.size()));
//                    advance(gs,ac);
//                    repairedChild.add(ac);
                }
            }
        }

        if(!(gs.getActiveTribeID() == getPlayerID())){return repairedChild;}

        while (!gs.isGameOver() && (gs.getActiveTribeID() == getPlayerID())){
            ArrayList<Action> allAvailableActions = this.allGoodActions(gs, m_rnd);
            Action a = allAvailableActions.get(m_rnd.nextInt(allAvailableActions.size()));
            advance(gs, a);
            repairedChild.add(a);
        }
        return repairedChild;
    }

    //give a random possible move as a mutation
    private Action mutation(GameState gs){

        ArrayList<Action> allAvailableActions = this.allGoodActions(gs, m_rnd);
        return  allAvailableActions.get(m_rnd.nextInt(allAvailableActions.size()));

    }

    public double eval(GameState gs, Individual actionSet){
        for(Action move : actionSet.getActions()){
            advance(gs, move);
        }
        actionSet.setGs(gs);
        return heuristic.evaluateState(gs);
    }

    private boolean checkActionFeasibility(Action a, GameState gs)
    {
        if(gs.isGameOver())
            return false;

        if(a instanceof UnitAction)
        {
            UnitAction ua = (UnitAction)a;
            int unitId = ua.getUnitId();
            Actor act = gs.getActor(unitId);
            if(!(act instanceof Unit) || act.getTribeId() != gs.getActiveTribeID())
                return false;
        }else if (a instanceof CityAction)
        {
            CityAction ca = (CityAction)a;
            int cityId = ca.getCityId();
            Actor act = gs.getActor(cityId);
            if(!(act instanceof City) || act.getTribeId() != gs.getActiveTribeID())
                return false;
        }

        boolean feasible = false;
        try{
            feasible = a.isFeasible(gs);
        }catch (Exception e) { }

        return feasible;
    }

    private void advance(GameState gs, Action move){
        this.fmCallsCount++;
        gs.advance(move,true);
    }
}
