import java.util.*;
import java.util.stream.Collectors;

class Cake {
    Cake previous; // предыдущие два куска, объединенные в один торт
    double k; // коэффициент масштабирования кусков
    double joint; // точка соединения двух кусков
    AllocCake alloc; // распределенная часть торта
    UnallocCake unalloc; // нераспределенная часть торта
    Cake() {
        k = 1;
        joint = 1;
        alloc = new AllocCake();
        unalloc = new UnallocCake();
    }
}
class AllocCake {
    double[][] prices; // стоимости частей 4х агентов по разным критериям
    double[] bonus; // стоимости бонуса для каждого из агентов по разным критериям
    List<int[]> maxIndex; // индексы
    int[] guaranteed; // индекс гарантированной части по стоимости для каждого агента
    Trim[] rm; // rightmost
    Trim[] srm; // secondRightmost
    AllocCake() {
        prices = new double[4][4];
        bonus = new double[4];
        maxIndex = new ArrayList<int[]>();
        for (int i = 0; i < 3; i++) maxIndex.add(new int[3]);
        guaranteed = new int[]{3,3,3,3};
        rm = new Trim[4];
        for (int i = 0; i < 4; i++) rm[i] = new Trim();
        srm = new Trim[4];
        for (int i = 0; i < 4; i++) srm[i] = new Trim();
    }
}
class UnallocCake {
    Section res1, res2; // два остатка после очередного дележа
    ArrayList<Integer> sp; // sp - significant piece
    Integer insp; // insp - insignificant piece
    UnallocCake() {
        res1 = new Section(0,1);
        res2 = null;
        sp = new ArrayList<>();
    }
}
class Section {
    double a, b; // начало и конец отрезка
    Section(double a, double b) {
        this.a = a;
        this.b = b;
    }
}
class Trim { // разрез
    double v; // значение
    int a; // агент, которому принадлежит разрез
    int iq; // индекс части, которой принадлежит разрез
}
enum Agents{A, B, C, D};

public class FairDivision {
    /*
    Представляем весь торт как отрезок [0,1]
    Разбиваем отрезок [0,1] на N неделимых кусочков.
    Каждый кусочек случайно оцениваем по четырем разным критериям.
     */
    static int N = 100;
    static double[][] values;
    static double[][] prices = new double[4][4];
    static List<Integer>[] dmnts = new List[4];

    public static void main(String[] args) {
        values = new double[4][N];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < N; j++) {
                values[i][j] = 1 - Math.random();
            }
        Cake cake = new Cake();
        for (int i = 0; i < 4; i++) {
            dmnts[i] = new ArrayList<>();
        }
        out:
        for (int cutter = 0; (cutter < 4) && (cake.unalloc.res1!=null); cutter++) {
            boolean dblDmntnFound = false;
            double[][] tableBonus = new double[4][4];
            ArrayList<Cake> cakes = new ArrayList<>();
            for (int i = 0; i < 4; i++) {
                cake = coreProtocol(cake, cutter);
                cakes.add(cake);
                if (cake.unalloc.res1==null) {
                    addPrices(cake.alloc.prices);
                    break out;
                }
                if (cake.previous != null) {
                    if (cake.alloc.rm[cake.unalloc.sp.get(0)].a != cake.previous.alloc.rm[cake.previous.unalloc.sp.get(0)].a) {
                        dmnts[cutter].add(cake.alloc.rm[cake.unalloc.sp.get(0)].a);
                        dmnts[cutter].add(cake.previous.alloc.rm[cake.previous.unalloc.sp.get(0)].a);
                        dblDmntnFound = true;
                        break;
                    }
                }
                for (int j = 0; j < 4; j++) {
                    tableBonus[i][j] = cake.alloc.bonus[j];
                }
            }
            if (!dblDmntnFound) {
                int minRow = 0;
                out1:
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        double sum = 0;
                        for (int k = 0; k < 4; k++) {
                            if (k != i) {
                                sum += tableBonus[k][j];
                            }
                        }
                        if (tableBonus[i][j] > sum) break;
                        if (j == 3) {
                            minRow = i;
                            break out1;
                        }
                    }
                }
                Cake mRCake = cakes.get(minRow);
                permutationProtocol(mRCake, cutter);
                dmnts[cutter].add(cake.alloc.rm[cake.unalloc.sp.get(0)].a);
                dmnts[cutter].add(mRCake.alloc.srm[mRCake.unalloc.sp.get(0)].a);
            }
            cake = coreProtocol(cake, cutter);
            cakes.add(cake);
            for (Cake c: cakes)
                addPrices(c.alloc.prices);
        }
        if (cake.unalloc.res1!=null) {
            postDoubleDominationProtocol(cake);
        }
        printPrices();
    }
    static void printPrices() {
        String[] agents = {"A","B","C","D"};
        System.out.printf("%20s%s\n"," ","Суммарная ценность частей, принадлежащих агенту");
        System.out.printf("%20s", " ");
        for (int i = 0; i < 4; i++) {
            System.out.printf("%-20s",agents[i]);
        }
        System.out.println();
        for (int i = 0; i < 4; i++) {
            System.out.printf("По оценке %s: ",agents[i]);
            for (int j = 0; j < 4; j++)
                System.out.printf("%-20.10f", prices[i][j]);
            System.out.println();
        }
    }
    static void addPrices(double[][] cakePrices) {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                prices[i][j] += cakePrices[i][j];
    }

    public static double cost(double x, double y, int agent, Cake cake) {
        if (x == y)
            return 0;
        double k = cake.k;
        double joint = cake.joint;
        if (k != 1) {
            double r1a = cake.previous.unalloc.res1.a;
            double r1b = cake.previous.unalloc.res1.b;
            double r2a = 0;
            if (cake.previous.unalloc.res2 != null ) {
                r2a = cake.previous.unalloc.res2.a;
            }
            if (x < joint && y > joint) {
                return cost(x * k + r1a, r1b, agent, cake.previous) +
                        cost(r2a, y * k + r2a - r1b + r1a, agent, cake.previous);
            } else if (y <= joint) {
                return cost(x * k + r1a,y * k + r1a, agent, cake.previous);
            } else {
                return cost(x * k + r2a - r1b + r1a,y * k + r2a - r1b + r1a, agent, cake.previous);
            }
        } else {
            double sum = 0;
            x *= N;
            y *= N;
            int begin = (int)Math.ceil(x);
            int end = (int)Math.floor(y);
            if (x != 0) {
                sum += values[agent][begin-1]*(begin-x);
            }
            for (int i = begin; i < end; i++)
                sum += values[agent][i];
            if ((int)y != N) {
                sum += values[agent][end]*(y-end);
            }
            return sum;
        }
    }

    public static double cut(double x, double cost, int agent, Cake cake) {
        if (cost == 0)
            return x;
        double k = cake.k;
        if (k != 1) {
            double r1a = cake.previous.unalloc.res1.a;
            double r1b = cake.previous.unalloc.res1.b;
            double r2a = cake.previous.unalloc.res2 != null ? cake.previous.unalloc.res2.a : 0;
            if (x < cake.joint) {
                return (cut(r1a + x*k, cost, agent, cake.previous) - r1a) / k;
            } else {
                double delta = x-cake.joint;
                return (cut(r2a+delta*k, cost, agent, cake.previous) - r2a + r1b-r1a) / k;
            }
        } else {
            double sum = 0;
            int end = (int)Math.ceil(x*N);
            if (x != 0) {
                if (values[agent][end-1] > cost) {
                    return (x*N + cost / values[agent][end-1]) / N;
                }
                sum += values[agent][end-1]*(end-x*N);
            }
            if (end >= cake.unalloc.res1.b*N) {
                for (; sum + values[agent][end] < cost; end++) {
                    sum += values[agent][end];
                }
                return (end + (cost - sum) / values[agent][end]) / N;
            } else {
                for (; sum + values[agent][end] < cost; end++) {
                    if (end >= cake.unalloc.res1.b*N) {
                        return cut(cake.unalloc.res2.a, cost-sum,agent,cake);
                    }
                    sum += values[agent][end];
                }
            }
            return (end + (cost - sum) / values[agent][end]) / N;
        }

    }

    private static Cake coreProtocol(Cake previousCake, int cutter) {
        Cake cake;
        Section res1 = previousCake.unalloc.res1;
        Section res2 = previousCake.unalloc.res2;
        if (res1.a != 0 || res1.b != 1) {
            cake = new Cake();
            cake.previous = previousCake;
            double deltaR1 = res1.b - res1.a;
            double deltaR2 = res2 != null ? res2.b - res2.a : 0;
            cake.k = deltaR1 + deltaR2;
            cake.joint = deltaR1 / cake.k;
        } else {
            cake = previousCake;
        }

        AllocCake alloc = cake.alloc;
        Trim[] rightmost = alloc.rm;
        Trim[] secondRightmost = alloc.srm;
        List<Integer> agents = new ArrayList<>(Arrays.asList(0,1,2,3));
        agents.remove((Integer)cutter);

        //cutter режет на 4 равные части
        double[] trims4 = new double[5];
        double priceQuarter = cost(0,1,cutter, cake)*0.25;
        trims4[0] = 0;
        for (int i = 1; i < 4; i++) {
            trims4[i] = cut(trims4[i-1], priceQuarter, cutter, cake);
        }
        trims4[4] = 1;

        //записываем стоимости 4-х частей для {1,2,3,4}
        double[][] pricesQuarters = new double[4][4];
        for (Integer i: agents) {
            for (int j = 0; j < 4; j++) {
                pricesQuarters[i][j] = cost(trims4[j], trims4[j+1], i, cake);
            }
        }
        for (int j = 0; j < 4; j++) {
            pricesQuarters[cutter][j] = priceQuarter;
        }

        //Находим индексы трех наиболее ценных частей для всех агентов кроме cutter'а
        List<int[]> maxIndex = alloc.maxIndex;
        List<Integer> mostQuarter = new ArrayList<>(Arrays.asList(0,1,2,3));
        for (int i = 0; i < 3; i++) {
            List<Double> list = Arrays.stream(pricesQuarters[agents.get(i)]).boxed().collect(Collectors.toList());
            maxIndex.set(i,list.stream().sorted(Collections.reverseOrder()).limit(3).mapToInt(list::indexOf).toArray());
            mostQuarter.remove((Integer)maxIndex.get(i)[0]);
        }

        //если все претендуют на разные доли, то отдаем каждому самую ценную часть
        if (mostQuarter.size() == 1) {
            for (int j = 0; j < 3; j++)
                for (int i = 0; i < 4; i++)
                    alloc.prices[i][agents.get(j)] = pricesQuarters[i][maxIndex.get(j)[0]];
            for (int i = 0; i < 4; i++)
                alloc.prices[i][cutter] = pricesQuarters[i][mostQuarter.get(0)];
            cake.unalloc.res1 = null;
            cake.unalloc.res2 = null;
            return cake;
        }

        //{I,J,K} приравнивают 1-ую и 2-ую по ценности части к 3-ей (каждый делает два разреза)
        List<Trim[]> trimsTo3 = new ArrayList<>();
        for (int i = 0; i < 3; i++) trimsTo3.add(new Trim[2]);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 2; j++)
                trimsTo3.get(i)[j] = new Trim();
        //количество разрезов в каждой части
        int[] numTrims = {0,0,0,0};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 2; j++) {
                trimsTo3.get(i)[j].v = cut(trims4[maxIndex.get(i)[j]],pricesQuarters[agents.get(i)][maxIndex.get(i)[j]] - pricesQuarters[agents.get(i)][maxIndex.get(i)[2]], agents.get(i), cake);
                trimsTo3.get(i)[j].iq = maxIndex.get(i)[j];
                trimsTo3.get(i)[j].a = agents.get(i);
                numTrims[maxIndex.get(i)[j]]++;
            }
        }

        //вычисляем 1ый и 2ой по величине правейшие разрезы каждой части
        for (int i = 0; i < 4; i++) {
            rightmost[i].v = trims4[i];
            rightmost[i].iq = i;
            rightmost[i].a = cutter;
            secondRightmost[i].v = trims4[i];
            secondRightmost[i].iq = i;
            secondRightmost[i].a = cutter;
        }

        ArrayList<Integer> agentsTo3 = new ArrayList<>(agents);
        ArrayList<Trim> rmTrims = new ArrayList<>();
        ArrayList<Integer> agentsTo2 = new ArrayList<>();
        for (int i = 0; i < 3; i++) {
            if (numTrims[maxIndex.get(i)[1]] == 1) {
                agentsTo2.add(agents.get(i));
            }
        }
        for (int i = 0; i < 3; i++)
            for (int j = i+1; j < 3; j++)
                if (maxIndex.get(i)[1] == maxIndex.get(j)[1] && numTrims[maxIndex.get(i)[1]] == 2 && agentsTo2.isEmpty()) {
                    agentsTo2.add(agents.get(i));
                    agentsTo2.add(agents.get(j));
                }

        //агенты, которые оценили одну и ту же часть как вторую по ценности
        //приравнивают 1-ую по ценности часть к 2-ой
        ArrayList<Trim> trimsTo2 = new ArrayList<>();
        for (int i = 0; i < agentsTo2.size(); i++) {
            trimsTo2.add(new Trim());
            trimsTo2.get(trimsTo2.size()-1).v = cut(trims4[maxIndex.get(i)[0]],pricesQuarters[agents.get(i)][maxIndex.get(i)[0]] - pricesQuarters[agents.get(i)][maxIndex.get(i)[1]], agents.get(i), cake);
            trimsTo2.get(trimsTo2.size()-1).iq = maxIndex.get(i)[0];
            trimsTo2.get(trimsTo2.size()-1).a = agents.get(i);
            cake.alloc.guaranteed[agents.get(i)] = 1;
        }
        //добавляем в общую выборку нужные разрезы
        for (int i = 0; i < agentsTo2.size(); i++) {
            rmTrims.add(trimsTo2.get(i));
        }
        agentsTo3.removeAll(agentsTo2);
        for (int i = 0; i < agentsTo3.size(); i++) {
            if (numTrims[maxIndex.get(agents.indexOf(agentsTo3.get(i)))[0]] != 1) {
                for (int j = 0; j < 2; j++) {
                    rmTrims.add(trimsTo3.get(agents.indexOf(agentsTo3.get(i)))[j]);
                }
            }
            else rmTrims.add(trimsTo3.get(agents.indexOf(agentsTo3.get(i)))[0]);
            cake.alloc.guaranteed[agentsTo3.get(i)] = 2;
        }
        for (Trim trim: rmTrims) {
            if (trim.v > rightmost[trim.iq].v) {
                secondRightmost[trim.iq] = rightmost[trim.iq];
                rightmost[trim.iq] = trim;
            } else if (trim.v > secondRightmost[trim.iq].v) {
                secondRightmost[trim.iq] = trim;
            }
        }
        for (int i = 0; i < 4; i++)
            for (int j = i+1; j < 4; j++) {
                if (rightmost[i].a == rightmost[j].a && rightmost[i].a != cutter) {
                    if (cost(secondRightmost[i].v,trims4[i+1], rightmost[i].a, cake)
                            > cost(secondRightmost[j].v,trims4[j+1], rightmost[i].a, cake)) {
                        rightmost[j] = secondRightmost[j];
                    } else {
                        rightmost[i] = secondRightmost[i];
                    }
                }
            }
        int rmCutter = 0;
        for (int i = 0; i < 4; i++)
            if (rightmost[i].a == cutter) rmCutter++;
        if (rmCutter > 1) {
            //находим оставшегося нережущего агента/-ов
            if (!agentsTo3.isEmpty()) {
                int sum = 0;
                for (int i = 0; i < 4; i++)
                    if (rightmost[i].a != cutter)
                        sum += rightmost[i].a;
                if (agentsTo2.contains(6 - cutter - sum)) {
                    rightmost[maxIndex.get(agents.indexOf(6 - cutter - sum))[1]].a = 6 - cutter - sum;
                } else {
                    rightmost[maxIndex.get(agents.indexOf(6 - cutter - sum))[2]].a = 6 - cutter - sum;
                }
            } else {
                for (int i = 0; i < 4; i++)
                    if (rightmost[i].a != cutter)
                        agentsTo2.remove((Integer)rightmost[i].a);
                for (int i = 0; i < agentsTo2.size(); i++) {
                    rightmost[maxIndex.get(i)[1]].a = agents.get(i);
                }
            }
        }
        //правейший агент каждой части получает кусок от правого края до 2го правейшего разреза
        Section[] res = {null, null};
        Integer[] resIndex = {null, null};
        for (int i = 0, k = 0; i < 4; i++) {
            if (secondRightmost[i].a != cutter) {
                for (int j = 0; j < 4; j++) {
                    alloc.prices[j][rightmost[i].a] = cost(secondRightmost[i].v,trims4[i+1], j, cake);
                }
                if (trims4[i] != secondRightmost[i].v) {
                    res[k] = new Section(trims4[i], secondRightmost[i].v);
                    resIndex[k] = i;
                    k++;
                }
            } else {
                for (int j = 0; j < 4; j++) {
                    alloc.prices[j][rightmost[i].a] = pricesQuarters[j][i];
                }
            }
            alloc.bonus[rightmost[i].a] = cost(secondRightmost[i].v,rightmost[i].v, rightmost[i].a, cake);
        }
        cake.unalloc.res1 = res[0];
        cake.unalloc.res2 = res[1];

        //определяем "значительную" часть
        double costR1 = cost(res[0].a,res[0].b,cutter,cake);
        double costR2 = 0;
        if (res[1] != null) {
            costR2= cost(res[1].a,res[1].b,cutter,cake);
        }
        if (costR1 > costR2) {
            cake.unalloc.sp.add(resIndex[0]);
            cake.unalloc.insp = resIndex[1];
        } else if (costR1 < costR2){
            cake.unalloc.sp.add(resIndex[1]);
            cake.unalloc.insp = resIndex[0];
        } else {
            cake.unalloc.sp.add(resIndex[0]);
            cake.unalloc.sp.add(resIndex[1]);
        }
        return cake;
    }

    static Cake permutationProtocol(Cake cake, int cutter) {
        AllocCake al = cake.alloc;
        Trim[] rm = al.rm;
        Trim[] srm = al.srm;
        Integer sp = cake.unalloc.sp.get(0);
        Integer insp = cake.unalloc.insp;
        List<Integer> agents = new ArrayList<>();
        for (int i = 0; i < 4; i++) agents.add(i);
        int I = rm[sp].a; // i := 1st
        int J = srm[sp].a; // j := 2nd
        agents.remove((Integer)I);
        agents.remove((Integer)J);
        agents.remove((Integer)cutter);
        int K = agents.get(0); // k := 3rd

        if (insp != null) {
            if (rm[insp].a == J) {
                if (srm[insp].a == I) {
                    swapCols(al.prices, I, J);
                } else if (srm[insp].a == J) {
                    swapCols(al.prices, I, J);
                } else {
                    swapCols(al.prices, I, J);
                    swapCols(al.prices, K, I);
                    if (al.prices[I][I] < al.prices[I][cutter]) {
                        swapCols(al.prices, I, cutter);
                    }
                }
            } else {
                if (rm[al.maxIndex.get(I)[al.guaranteed[I]]].a == J) {
                    swapCols(al.prices, I, J);
                } else if (rm[al.maxIndex.get(I)[al.guaranteed[I]]].a == cutter) {
                    swapCols(al.prices, I, J);
                    swapCols(al.prices, I, cutter);
                } else {
                    swapCols(al.prices, I, J);
                    swapCols(al.prices, I, K);
                    if (al.prices[K][K] < al.prices[I][cutter]) {
                        swapCols(al.prices, K, cutter);
                    }
                }
            }
        } else {
            if (rm[al.maxIndex.get(I)[al.guaranteed[I]]].a == J) {
                swapCols(al.prices, I, J);
            } else if (rm[al.maxIndex.get(I)[al.guaranteed[I]]].a == cutter) {
                swapCols(al.prices, I, J);
                swapCols(al.prices, I, cutter);
            } else {
                swapCols(al.prices, I, J);
                swapCols(al.prices, I, K);
                if (al.prices[K][K] < al.prices[I][cutter]) {
                    swapCols(al.prices, K, cutter);
                }
            }
        }
        return cake;
    }

    static void swapCols(double[][] m, int k, int l) {
        for (int i = 0; i < m.length; i++) {
            double t = m[i][k];
            m[i][k] = m[i][l];
            m[i][l] = t;
        }
    }

    static void postDoubleDominationProtocol(Cake cake) {
        Section res1 = cake.unalloc.res1;
        Section res2 = cake.unalloc.res2;

        Cake allResidue = new Cake();
        allResidue.previous = cake;
        double deltaR1 = res1.b - res1.a;
        double deltaR2 = res2 != null ? res2.b - res2.a : 0;
        allResidue.k = deltaR1 + deltaR2;
        allResidue.joint = deltaR1 / allResidue.k;

        List<Integer> agents = new ArrayList<>(Arrays.asList(1,2,3));
        Integer I = 0; //1st
        Integer K = dmnts[I].get(0); //3rd
        Integer M = dmnts[I].get(1); //4th
        agents.remove(K);
        agents.remove(M);
        Integer J = agents.get(0); //2nd

        if (dmnts[J].containsAll(dmnts[I])) {
            //K and M divide the residue by divideAndChoose
            divideAndChoose(allResidue, K, M);
        } else {
            if (dmnts[J].contains(K)) {
                if (dmnts[M].contains(K)) {
                    //give all to K
                    for (int i = 0; i < 4; i++) {
                        prices[i][K] += cost(0, 1, i, allResidue);
                    }
                } else {
                    //K cuts then pieces are given in order I,J,M
                    divideByOrder(allResidue, K, I, J, M);
                }
            } else {
                if (dmnts[K].contains(M)) {
                    //give all to M
                    for (int i = 0; i < 4; i++) {
                        prices[i][M] += cost(0, 1, i, allResidue);
                    }
                } else {
                    //M cuts then pieces are given in order I,J,K
                    divideByOrder(allResidue, M, I, J, K);
                }
            }
        }
    }

    static void divideAndChoose(Cake cake, int I, int J) {
        double priceHalf = cost(0,1, I, cake) * 0.5;
        double trim = cut(0, priceHalf, I, cake);
        if (cost(0, trim, J, cake) > cost(trim, 1, J, cake)) {
            for (int i = 0; i < 4; i++) {
                prices[i][J] += cost(0, trim, i, cake);
                prices[i][I] += cost(trim, 1, i, cake);
            }
        } else {
            for (int i = 0; i < 4; i++) {
                prices[i][I] += cost(0, trim, i, cake);
                prices[i][J] += cost(trim, 1, i, cake);
            }
        }
    }

    static void divideByOrder(Cake cake, int cutter, int I, int J, int K) {
        double[] trims4 = new double[5];
        double priceQuarter = cost(0,1,cutter, cake)*0.25;
        trims4[0] = 0;
        for (int i = 1; i < 4; i++) {
            trims4[i] = cut(trims4[i-1], priceQuarter, cutter, cake);
        }
        trims4[4] = 1;

        List<Integer> agents = new ArrayList<>(Arrays.asList(I,J,K));
        double[][] pricesQuarters = new double[4][4];
        for (Integer i: agents) {
            for (int j = 0; j < 4; j++) {
                pricesQuarters[i][j] = cost(trims4[j], trims4[j+1], i, cake);
            }
        }
        for (int j = 0; j < 4; j++) {
            pricesQuarters[cutter][j] = priceQuarter;
        }

        List<Integer> quarterIndex = new ArrayList<>(Arrays.asList(0,1,2,3));
        for (Integer agent : agents) {
            double max = 0;
            int maxj = 0;
            for (Integer index : quarterIndex) {
                if (pricesQuarters[agent][index] > max) {
                    max = pricesQuarters[agent][index];
                    maxj = index;
                }
            }
            quarterIndex.remove((Integer) maxj);
            for (int k = 0; k < 4; k++) {
                prices[k][agent] += pricesQuarters[k][maxj];
            }
        }
        for (int k = 0; k < 4; k++) {
            prices[k][cutter] += pricesQuarters[k][quarterIndex.get(0)];
        }
    }
}
