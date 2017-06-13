package com.company;

import java.util.ArrayList;

/**
 * Created by ethur on 5/12/17.
 */
 class Gene {

    // transcript
    private String chromosome;
    private int start;
    private int stop;
    private char orientation;
    private String name;
    private double score;

    private ArrayList<Position> positionList = new ArrayList<>();

    private String description;

    //chromosome, name, start, stop, score ,orientation, description


     Gene(String chromosome, int start, int stop, char orientation, String name, double score) {
        this.chromosome = chromosome;
        this.start = start;
        this.stop = stop;
        this.orientation = orientation;
        this.name = name;
        this.score = score;
    }

    public ArrayList<Position> getPositionList() {
        return positionList;
    }

     String getChromosome() {
        return chromosome;
    }

     int getStart() {
        return start;
    }

     int getStop() {
        return stop;
    }

     char getOrientation() {
        return orientation;
    }

    String getName() {
        return name;
    }

    public double getScore() {
        return score;
    }


    void addToPositionList(Position pos) {
        positionList.add(pos);
    }


    int qualityScoring (){

        int totalQuality = 0;
        for (Position pos: positionList
             ) {
            totalQuality += pos.getCoverage();
        }
        //System.out.println("[STATUS] Gene Quality total " + totalQuality);
        return totalQuality;

    }


    // ToDo  orientation specificity
    boolean containsPosition(Position position) {

        if (chromosome.equals(position.getChromosome()) && orientation == position.getOrientation()) {
            System.out.println("right area");
            if (start < position.getPosition() && stop > position.getPosition()) {
                return true;
            }

        }


        return false;
    }


    @Override
    public String toString() {

        int[] arrayCoverage = new int[(stop - start)];  // array of length gene


        String resultString = "";

        StringBuilder stringBldr = new StringBuilder((stop - start));

        for (Position pos : positionList
                ) {
            int relativePosition = pos.getPosition() - (start+1);
            arrayCoverage[relativePosition] = pos.getCoverage();
        }

        for (int coverage : arrayCoverage
                ) {
            stringBldr.append(coverage);
            stringBldr.append(";");
        }

        return name + '\t' + stringBldr.toString();
    }
}
