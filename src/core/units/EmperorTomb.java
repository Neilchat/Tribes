package core.units;

import core.Types;
import core.game.Board;

public class EmperorTomb extends Building {

    private final Types.BUILDING TYPE = Types.BUILDING.EMPERORS_TOMB;
    private final int PRODUCTION = 3;
    private final int points = 400;
    private final Types.TERRAIN TERRAIN_CONSTRAINT_1 = Types.TERRAIN.PLAIN;
    private final Types.TERRAIN TERRAIN_CONSTRAINT_2 = Types.TERRAIN.SHALLOW_WATER;

    public EmperorTomb(int x, int y) {
        super(x, y);
    }

    @Override
    public Building copy() {
        return new EmperorTomb(getX(), getY());
    }

    @Override
    public Types.BUILDING getTYPE() {
        return TYPE;
    }

    @Override
    public int getCOST() {
        return 0;
    }

    @Override
    public int getPRODUCTION() {
        return PRODUCTION;
    }

    @Override
    public int getPoints() {
        return points;
    }

    @Override
    public boolean is_buildable(Board board) {
        return board.getTerrainAt(getX(), getY()).equals(TERRAIN_CONSTRAINT_1) || board.getTerrainAt(getX(), getY()).equals(TERRAIN_CONSTRAINT_2);
    }
}
