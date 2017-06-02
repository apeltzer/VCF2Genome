package VCFTools;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by peltzer on 15/01/2017.
 * Calculated ambigous characters
 *
 */
public class AmbigoutyCalculator {

    public AmbigoutyCalculator() {
    }


    public char getAmbiguousBase(char c1, char c2) {
        Set<Character> chars = new HashSet<Character>();
        chars.add(c1);
        chars.add(c2);

        if (c1 == c2)
            return c1;

        if (chars.contains('G') && chars.contains('A'))
            return 'R';
        if (chars.contains('T') && chars.contains('C'))
            return 'Y';
        if (chars.contains('G') && chars.contains('T'))
            return 'K';
        if (chars.contains('A') && chars.contains('C'))
            return 'M';
        if (chars.contains('G') && chars.contains('C'))
            return 'S';
        if (chars.contains('A') && chars.contains('T'))
            return 'W';

        System.err.println("Illegal arguments in function getAmbiguousBase: " + c1 + ", " + c2);

        return 'N';
    }
}
