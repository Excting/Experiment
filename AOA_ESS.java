

import java.io.File;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.Date;
import jxl.Workbook;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;



public class AOA_ESS {
	private static final int RUN = 30;
	private static final int pop_num = 50;

	private static final double miu = 0.499;
	private static final int a = 5;
	private double MOP_Max=1.0;
	private double MOP_Min=0.2;
	private double eps=0.00000000000001;

	private static final int K = 10;
	private static final double alpha = 0.001;


	private static final double Pc = 0.2;
	private static final double F = 0.5;


	private static int fun_num;

	private static int R;

	private static int PATHNUM;

	private static int NODENUM;

	private static int[][]infection;
	private static final int BRANCH = 15;

	private static boolean[][] visit;

	private static int[] record = new int[pop_num];

	private static String[] PATH;

	private static final int step_length = 2;


	static double start;
	static double finish;
	static double[] runtime;
	static double[] coverage;
	static int[] case_num;
	static int[] obj = new int[1];


	private static int[] lb;
	private static int[] ub;

	private static final int MCN = 300000;
	private static int col;

	public static void main(String[] args){

		for(fun_num = 7;fun_num< 13;fun_num++)
		{
			System.out.println("FUNCTION = " + fun_num);

			setFunctionParameters();

			lb = new int[R];
			ub = new int[R];
			visit = new boolean[NODENUM][BRANCH];
			PATH = new String[PATHNUM];
			runtime = new double[RUN];
			coverage = new double[RUN];
			case_num = new int[RUN];
			infection = new int[R][NODENUM];
			setDomainAndEncoding();

			for (int run = 0; run < RUN; run++)
			{

				for (int i = 0; i < NODENUM; i++)
					for (int j = 0; j < BRANCH; j++)
						visit[i][j] = false;


				int[][] x = new int[pop_num][R];
				int[][] x_re = new int[pop_num][R];
				int[][] v = new int[pop_num][R];
				int[][] v_new = new int[pop_num][R];

				double[] fitness_x = new double[pop_num];
				double[] fitness_x_re = new double[pop_num];
				double[] fitness_v = new double[pop_num];
				double[] fitness_v_new = new double[pop_num];
				boolean[] status = new boolean[PATHNUM];

				int[] res = new int[PATHNUM];
				int[][] solution = new int[PATHNUM][R];

				int path;
				obj[0] = 0;

				Date mydate = new Date();
				start = mydate.getTime();


				for (int i = 0; i < pop_num; i++) {
					for (int j = 0; j < R; j++) {
						x[i][j] = (int) (lb[j] + Math.random() * (ub[j] - lb[j]));
					}

					if (obj[0] == PATHNUM)
						break;

					path = pathnum(x[i], fun_num);
					record[i] = path;
					if (!status[path]) {
						for (int j = 0; j < R; j++)
							solution[path][j] = x[i][j];
						status[path] = true;
						obj[0]++;
						res[path] = 0;
						nodeiscoverage(x[i], fun_num);
					}
					case_num[run] = case_num[run] + 1;
				}


				for(int i = 0; i < pop_num; i++){
					for(int j = 0; j < R; j++){
						x_re[i][j] = (int)(lb[j]+ub[j])-x[i][j];
						if(x_re[i][j]<lb[j]||x_re[i][j]>=ub[j]){
							x_re[i][j] = (int) (lb[j] + Math.random() * (ub[j] - lb[j]));
						}
					}
					if (obj[0] == PATHNUM)
						break;
					path = pathnum(x_re[i], fun_num);
					record[i] = path;

					if (!status[path]) {
						for (int j = 0; j < R; j++)
							solution[path][j] = x_re[i][j];
						status[path] = true;
						obj[0]++;
						nodeiscoverage(x_re[i], fun_num);
						res[path] = 0;
					}
					case_num[run] = case_num[run] + 1;
				}


				for (int i = 0; i < pop_num; i++) {

					fitness_x[i] = benchmarkfunction(x[i], fun_num, -1);
					fitness_x_re[i] = benchmarkfunction(x_re[i], fun_num, -1);
					case_num[run] = case_num[run] + 2;

					if (fitness_x_re[i] > fitness_x[i])
					{
						for (int j = 0; j < R; j++)
							x[i][j] = x_re[i][j];
						fitness_x[i] = fitness_x_re[i];
					}
					if (obj[0] == PATHNUM) {
						break;
					}
				}

				while (case_num[run] <= MCN && obj[0] < PATHNUM)
				{
					AOA_search(x,v_new,fitness_x, fitness_v_new,solution, status, run, obj,case_num[run],res);
					ScatterSearch(x, fitness_x, solution, status, run, obj,res);
					if(obj[0] == PATHNUM)
						break;
				}

				Date mydate2 = new Date();
				finish = mydate2.getTime();
				runtime[run] = finish - start;
				System.out.println();
				System.out.println("运行时间=" + runtime[run] + "ms");

				coverage[run] = obj[0] * 100.0 / PATHNUM;
				System.out.println("路径覆盖率=" + coverage[run] + "%");
				System.out.println("最优解为：");

				for (int k = 0; k < PATHNUM; k++)
				{
					if (status[k]) {
						System.out.print("path" + k + ":");
						for (int j = 0; j < R; j++)
							System.out.print((int) solution[k][j] + " ");
						System.out.println();
					} else
						System.out.println("path" + k + "没被覆盖.");
				}
				Arrays.sort(res);
				System.out.print("[");
				for(int kk = 0; kk<PATHNUM;kk++){
					System.out.print((int)res[kk]+",");
				}
				System.out.print("]");
				System.out.println();

			}

			double time_sum = 0, time_average, coverage_sum = 0, coverage_average, case_average;
			int case_sum = 0;
			for (int run = 0; run < RUN; run++) {
				time_sum = time_sum + runtime[run];
				coverage_sum = coverage_sum + coverage[run];
				case_sum = case_sum + case_num[run];
			}
			time_average = time_sum / RUN;
			coverage_average = coverage_sum / RUN;
			case_average = case_sum / RUN;

			System.out.println("time_sum = " + time_sum + "ms");
			System.out.println("time_average = " + time_average + "ms");
			System.out.println("case_sum = " + case_sum);
			System.out.println("case_average = " + case_average);
			System.out.println("coverage_sum = " + coverage_sum + "%");
			System.out.println("coverage_average = " + coverage_average + "%");

			try
			{
				File file = new java.io.File("F:/paper/data/AOA_ESS.xls");

				Workbook book = Workbook.getWorkbook(file);
				WritableWorkbook wbook = Workbook.createWorkbook(file, book);
				WritableSheet sheet = wbook.getSheet(0);

				jxl.write.Number ID = new jxl.write.Number(col, 0, fun_num);
				sheet.addCell(ID);

				for (int run = 0; run < RUN; run++) {
					int q = run;
					jxl.write.Number number = new jxl.write.Number(col, q+1,case_num[run]);
					jxl.write.Number number2 = new jxl.write.Number(col, q+RUN+10,coverage[run]);
					jxl.write.Number number3 = new jxl.write.Number(col, q+15+RUN*2, runtime[run]);
					sheet.addCell(number);
					sheet.addCell(number2);
					sheet.addCell(number3);
				}

				double case_ave = getAverage(case_num, RUN);
				jxl.write.Number number1 = new jxl.write.Number(col, RUN + 4, case_ave);
				sheet.addCell(number1);
				double case_std = getStandardDevition(case_num, RUN);
				jxl.write.Number number2 = new jxl.write.Number(col, RUN + 5, case_std);
				sheet.addCell(number2);
				jxl.write.Number number3 = new jxl.write.Number(col, RUN + 6, time_sum);
				sheet.addCell(number3);

				wbook.write();
				wbook.close();

			} catch (Exception e) {
				System.out.println(e);
			}
		}
	}

	private static int jugde(String s, String s1) {
		int cut = 0;
		for(int i = 0; i  <s.length();i++){
			if(s.charAt(i)==s1.charAt(i)){
				cut++;
			}else{
				return cut;
			}
		}
		return  cut;
	}

	public static void AOA_search(int[][] x,int[][] v_new, double[] fitness_x,double[] fitness_v_new, int[][] solution, boolean[] status, int run, int[] obj, int citer,int[] res) {
		double MOP_Max=1.0;
		double MOP_Min=0.2;
		double eps=0.00000000000001;
		double MOP=1.0-(Math.pow((double)citer,1/a)/Math.pow(MCN,1/a));
		double MOA=MOP_Min+(double)citer*((MOP_Max-MOP_Min)/(double)MCN);
		for (int i = 0; i < pop_num; i++)
		{

			int target_path = random_UncoverPath(status, record[i]);
			int path;
			for(int j = 0; j < R;j++){

				int best;
				double[] fitness_temp = new double[step_length+1];

				best = x[i][j];
				int firstPath = pathnum(x[i], fun_num);
				if(Math.random()<MOA){
					if(Math.random()<0.5){
						v_new[i][j] = (int) (best/(MOP+eps)*((ub[j]-lb[j])*miu+lb[j]));
					}else{
						v_new[i][j] = (int) (best*(MOP+eps)*((ub[j]-lb[j])*miu+lb[j]));
					}
				}else{
					if(Math.random()<0.5){
						v_new[i][j] = (int) (best-(MOP)*((ub[j]-lb[j])*miu+lb[j]));
					}else{
						v_new[i][j] = (int) (best+(MOP)*((ub[j]-lb[j])*miu+lb[j]));
					}
				}
			}
			path = pathnum(v_new[i], fun_num);
			record[i] = path;
			if (!status[path]) {
				for (int j = 0; j < R; j++)
					solution[path][j] = v_new[i][j];
				status[path] = true;
				obj[0]++;
				res[path] = case_num[run];
				nodeiscoverage(v_new[i], fun_num);
			}

			if (obj[0] == PATHNUM)
				break;
		}
		for (int i = 0; i < pop_num; i++) {

			fitness_x[i] = benchmarkfunction(x[i], fun_num, -1);
			fitness_v_new[i] = benchmarkfunction(v_new[i], fun_num, -1);
			case_num[run] = case_num[run] + 1;

			if (fitness_v_new[i] > fitness_x[i])
			{
				for (int j = 0; j < R; j++)
					x[i][j] = v_new[i][j];
				fitness_x[i] = fitness_v_new[i];
			}
			if (obj[0] == PATHNUM) {
				break;
			}
		}

	}

public static void ScatterSearch(int[][]x, double[]fitness_x, int[][] solution, boolean[] status, int run, int[] obj,int[] res)
{
	for(int i=0;i<pop_num;i++)
	{

		int target_path = random_UncoverPath(status,record[i]);
		int path;


		for(int dim=0;dim<R;dim++)
		{

			int j = dim;
			int best;
			double[] fitness_temp = new double[step_length+1];

			int step;
			best = x[i][j];

			int firstPath = pathnum(x[i], fun_num);

			step = (ub[j]-lb[j])/step_length;

			while(step>1){

				int[] temp = getIndex(lb[j],ub[j],best,step);
				int[] temp_re  = new int[temp.length];
				for(int te = 0; te < temp.length; te++){
						 if(temp[te]<lb[j]||temp[te]>ub[j]){
							 temp_re[te] = lb[j]+(int)(Math.random()*(ub[j]-lb[j]));
						 }else{
							 temp_re[te] = lb[j]+ub[j] - temp[te];
						 }

				}
				for(int k=0;k<step_length;k++){

					x[i][j] = temp[k];
					if(x[i][j]>ub[j]||x[i][j]<lb[j])
					{
							continue;
					}


					path = pathnum(x[i], fun_num);

					update_Infection(firstPath, path, j);

					if(!status[path])
					{
						for(int t=0;t<R;t++)
							solution[path][t] = x[i][t];
						status[path] = true;
						obj[0]++;
						nodeiscoverage(x[i],fun_num);
					}
					if(obj[0] == PATHNUM)
						break;

					if(status[target_path])
						target_path = random_UncoverPath(status,path);


					int[] z = x[i];
					z[j] = temp_re[k];

					path = pathnum(z, fun_num);

					update_Infection(firstPath, path, j);

					if(!status[path])
					{
						for(int t=0;t<R;t++)
							solution[path][t] = z[t];
						status[path] = true;
						obj[0]++;
						nodeiscoverage(z,fun_num);
					}



					if(obj[0] == PATHNUM)
						break;

					if(status[target_path])
						target_path = random_UncoverPath(status,path);

					double a = benchmarkfunction(x[i], fun_num, target_path);
					double b = benchmarkfunction(z, fun_num, target_path);
					fitness_temp[k] = Math.max(a,b);
					if(b>a){
						x[i] = z;
					}
					case_num[run] = case_num[run] + 1;
				}
				int best_index = getBestIndex(fitness_temp);
				x[i][j] = temp[best_index];
				best = temp[best_index];
				fitness_x[i] = fitness_temp[best_index];

				step = step/step_length;
			}

			step = 1;
			int[] temp = getIndex(lb[j],ub[j],best,step);

			for(int k=0;k<step_length+1;k++)
			{
				x[i][j] = temp[k];
				if(x[i][j]>ub[j]||x[i][j]<lb[j])
					continue;

				path = pathnum(x[i], fun_num); //
				update_Infection(firstPath, path, j);

				if (!status[path]) {
					for (int t = 0; t < R; t++)
						solution[path][t] = x[i][t];
					status[path] = true;
					obj[0]++;
					nodeiscoverage(x[i], fun_num);
					res[path] = case_num[run];

				}

				if(obj[0] == PATHNUM)
					break;

				if(status[target_path])
					target_path = random_UncoverPath(status,path);



				if(status[target_path])
					target_path = random_UncoverPath(status,path);

				fitness_temp[k]  = benchmarkfunction(x[i], fun_num, target_path);

				case_num[run] = case_num[run] + 1;
			}

			int best_index = getBestIndex(fitness_temp);
			x[i][j] = temp[best_index];
			fitness_x[i] = fitness_temp[best_index];
		}
		if(obj[0] == PATHNUM)
			break;
	}
}
	public static void update_Infection(int firstPath, int path, int j)
	{
		if(firstPath!=path){
			for(int length=0;length<NODENUM;length++)
				if(PATH[firstPath].charAt(length)!=PATH[path].charAt(length)){
					if(PATH[firstPath].charAt(length)!=' '&&PATH[path].charAt(length)!=' ')
						infection[j][length]++;
				}
		}

	}

	public static int random_UncoverPath(boolean[] status,int path) {

		int[] similar = new int[PATHNUM];
		for(int i=0;i<PATHNUM;i++)
		{
			for(int j=0;j<NODENUM;j++)
			{
				if(PATH[path].charAt(j)==PATH[i].charAt(j))
					similar[i]++;
			}
		}
		int[] possible = new int[PATHNUM];
		for(int i=0;i<PATHNUM;i++)
			if(!status[i])
				possible[i] = similar[i];
			else possible[i] = 0;

		int temp=0; int index=0;
		for(int i=0;i<PATHNUM;i++)
			temp+=possible[i];	//统计总数
		if (temp == 0) {
			index = (int)(Math.random()*PATHNUM);
			while(status[index]){
				index = (int)(Math.random()*PATHNUM);
			}
		}else{

		}
		return index;
	}

	public static int checksum(String ISBN){
		int sum=0; int k=0;
		for(int i=0;i<ISBN.length();i++){
			if(ISBN.charAt(i)-'0'>=0&&ISBN.charAt(i)-'0'<=9){
				k++;
				sum+=(ISBN.charAt(i)-'0')*k;
			}
			if(ISBN.charAt(i)=='X'||ISBN.charAt(i)=='x')
			{
				k++;
				sum+=10*k;
			}
		}
		return sum;
	}


	public static void setFunctionParameters(){
		if (fun_num == 1) {
			R = 3;
			PATHNUM = 2;
			NODENUM = 1;
			col = 1;
		}
		if (fun_num == 2) {
			R = 2;
			PATHNUM = 9;
			NODENUM = 5;
			col = 2;
		}
		if (fun_num == 3) {
			R = 7;
			PATHNUM = 9;
			NODENUM = 6;
			col = 3;
		}
		if (fun_num == 4) {
			R = 7;
			PATHNUM = 5;
			NODENUM = 3;
			col = 4;
		}
		if (fun_num == 5) {
			R = 5;
			PATHNUM = 6;
			NODENUM = 3;
			col = 5;
		}
		//
		//static int[][] G={{1,3,2,1},{2,2,9,5},{3,7,9,6},{4,7,5,3},{5,5,6,3},{6,8,7,4},{7,7,48,4},{8,6,3,2},{9,11,12,7},{10,4,3,2},{11,8,4,2},{12,6,10,5}};
		if (fun_num == 6) {
			R = 8;
			PATHNUM = 7;
			NODENUM = 4;
			col = 6;
		}

		if (fun_num == 7) {
			R = 7;
			PATHNUM = 48;
			NODENUM = 4;
			col = 13;
		}
		if (fun_num == 8) {
			R = 6;
			PATHNUM = 3;
			NODENUM = 2;
			col = 14;
		}
		if (fun_num == 9) {
			R = 11;
			PATHNUM = 12;
			NODENUM = 7;
			col = 15;
		}
		if (fun_num == 10) {
			R = 4;
			PATHNUM = 3;
			NODENUM = 2;
			col = 16;
		}
		if (fun_num == 11) {
			R = 8;
			PATHNUM = 4;
			NODENUM = 2;
			col = 17;
		}
		if (fun_num == 12) {
			R = 6;
			PATHNUM = 10;
			NODENUM = 5;
			col = 18;
		}
	}

	public static void setDomainAndEncoding() {
		// transmit(1,3,2,1) 100%;
		if (fun_num == 1) {
			for (int i = 0; i < R; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = 255;
			}
		} // send(2,2,9,5)66%;
		if (fun_num == 2) {
			for (int i = 0; i < R; i++) {
				lb[i] = -1000000000; // 初始化上下界
				ub[i] = 1000000000;
			}
		} // processEvent(3,7,9,6)100%;
		if (fun_num == 3) {
			lb[0] = 0;
			ub[0] = 100;
			lb[1] = 0;
			ub[1] = 4;
			lb[2] = -2;
			ub[2] = 10000;
			lb[3] = 0;
			ub[3] = 2;
			lb[4] = 0;
			ub[4] = 10000;
			lb[5] = 0;
			ub[5] = 10000;
			lb[6] = -2;
			ub[6] = 100000;
		} // executeTuple(4,7,5,3)100%;
		if (fun_num == 4) {
			lb[0] = 1;
			ub[0] = 3;
			for (int i = 1; i < R; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = 255;
			}
		} // checkCloudletCompletion(5,5,6,3)100%;
		if (fun_num == 5) {
			lb[0] = 0;
			ub[0] = 2;
			for (int i = 1; i < R - 1; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = 255;
			}
			lb[R - 1] = 0;
			ub[R - 1] = 2;
		} // getResultantTuple(6,8,7,4)100%*/
		if (fun_num == 6) {
			for (int i = 0; i < 6; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = 255;
			}
			lb[6] = 0;
			ub[6] = 2;
			lb[7] = 1;
			lb[7] = 4;
		}
		// initFactory(13,7,48,4)
		if (fun_num == 7) {
			for (int i = 0; i < 6; i++) {
				lb[i] = 0;
				ub[i] = Integer.MAX_VALUE;
			}
			lb[6] = 0;
			ub[6] = Integer.MAX_VALUE;
			// ub[6] = 100;
		} // CleanXmlAnnotator(14,6,3,2)
		if (fun_num == 8) {
			for (int i = 0; i < R; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = Integer.MAX_VALUE;
			}
		} // WordsToSentencesAnnotator(15,11,12,7)
		if (fun_num == 9) {
			lb[0] = 0;
			ub[0] = 1;
			lb[1] = 0;
			ub[1] = 1;
			for (int i = 2; i < R; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = Integer.MAX_VALUE;
			}
		} // annotate(16,4,3,2)
		if (fun_num == 10) {
			for (int i = 0; i < R - 1; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = Integer.MAX_VALUE;
			}
			lb[R - 1] = 1;
			ub[R - 1] = Integer.MAX_VALUE;
		} // NERClassifierCombiner(17,8,4,2)
		if (fun_num == 11) {
			for (int i = 0; i < R - 1; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = Integer.MAX_VALUE;
			}
			lb[R - 1] = 0;
			ub[R - 1] = 1;
		} // setTrueCaseText(18,6,10,5)
		if (fun_num == 12) {
			for (int i = 0; i < R - 1; i++) {
				lb[i] = 0; // 初始化上下界
				ub[i] = Integer.MAX_VALUE;
			}
			lb[R - 1] = 0;
			ub[R - 1] = 1;
		} // Triangle(3,4,4);Factorial(1,2,1);BubbleSort(10,2,1);GCD(2,4,2);Middle(3,4,3);Commission(3,3,2);decision(4,11,4)


		switch (fun_num) {
			case 1:
				PATH[0] = "0";
				PATH[1] = "1";
				break;
			case 2:
				PATH[0] = "0    ";
				PATH[1] = "100  ";
				PATH[2] = "1010 ";
				PATH[3] = "10110";
				PATH[4] = "10111";
				PATH[5] = "110  ";
				PATH[6] = "1110 ";
				PATH[7] = "11110";
				PATH[8] = "11111";
				break;
			case 3:
				PATH[0] = "0     ";
				PATH[1] = "10    ";
				PATH[2] = "110   ";
				PATH[3] = "111 00";
				PATH[4] = "111 01";
				PATH[5] = "111 1 ";
				PATH[6] = "12 0  ";
				PATH[7] = "12 1  ";
				PATH[8] = "13    ";
				break;
			case 4:
				PATH[0] = "000";
				PATH[1] = "001";
				PATH[2] = "010";
				PATH[3] = "011";
				PATH[4] = "1  ";
				break;
			case 5:
				PATH[0] = "000";
				PATH[1] = "001";
				PATH[2] = "010";
				PATH[3] = "011";
				PATH[4] = "1 0";
				PATH[5] = "1 1";
				break;
			case 6:
				PATH[0] = "0000";
				PATH[1] = "0001";
				PATH[2] = "001 ";
				PATH[3] = "0100";
				PATH[4] = "0101";
				PATH[5] = "011 ";
				PATH[6] = "1   ";
				break;
			case 7:
				PATH[0] = "0000";
				PATH[1] = "0001";
				PATH[2] = "0002";
				PATH[3] = "0003";
				PATH[4] = "0004";
				PATH[5] = "0005";
				PATH[6] = "0006";
				PATH[7] = "0007";

				PATH[8] = "0010";
				PATH[9] = "0011";
				PATH[10] = "0012";
				PATH[11] = "0013";
				PATH[12] = "0014";
				PATH[13] = "0015";
				PATH[14] = "0016";
				PATH[15] = "0017";

				PATH[16] = "01 0";
				PATH[17] = "01 1";
				PATH[18] = "01 2";
				PATH[19] = "01 3";
				PATH[20] = "01 4";
				PATH[21] = "01 5";
				PATH[22] = "01 6";
				PATH[23] = "01 7";

				PATH[24] = "1000";
				PATH[25] = "1001";
				PATH[26] = "1002";
				PATH[27] = "1003";
				PATH[28] = "1004";
				PATH[29] = "1005";
				PATH[30] = "1006";
				PATH[31] = "1007";

				PATH[32] = "1010";
				PATH[33] = "1011";
				PATH[34] = "1012";
				PATH[35] = "1013";
				PATH[36] = "1014";
				PATH[37] = "1015";
				PATH[38] = "1016";
				PATH[39] = "1017";

				PATH[40] = "11 0";
				PATH[41] = "11 1";
				PATH[42] = "11 2";
				PATH[43] = "11 3";
				PATH[44] = "11 4";
				PATH[45] = "11 5";
				PATH[46] = "11 6";
				PATH[47] = "11 7";
				break;
			case 8:
				PATH[0] = "00";
				PATH[1] = "01";
				PATH[2] = "1 ";
				break;
			case 9:
				PATH[0] = "000    ";
				PATH[1] = "001    ";
				PATH[2] = "01     ";
				PATH[3] = "1  0   ";
				PATH[4] = "1  1000";
				PATH[5] = "1  1001";
				PATH[6] = "1  1010";
				PATH[7] = "1  1011";
				PATH[8] = "1  1100";
				PATH[9] = "1  1101";
				PATH[10] = "1  1110";
				PATH[11] = "1  1111";
				break;
			case 10:
				PATH[0] = "00";
				PATH[1] = "01";
				PATH[2] = "1 ";
				break;
			case 11:
				PATH[0] = "00";
				PATH[1] = "01";
				PATH[2] = "10";
				PATH[3] = "11";
				break;
			case 12:
				PATH[0] = "0   0";
				PATH[1] = "0   1";
				PATH[2] = "10  0";
				PATH[3] = "10  1";
				PATH[4] = "110 0";
				PATH[5] = "110 1";
				PATH[6] = "11100";
				PATH[7] = "11101";
				PATH[8] = "11110";
				PATH[9] = "11111";
				break;
		}
	}

	public static void nodeiscoverage(int[] x , int func_num)
	{
		if(func_num == 1)
		{
			char tuple[]= new char[R];
			for(int i=0;i<R;i++)
				tuple[i] = (char) x[i];

			if(tuple[0]=='O'&&tuple[1]=='l'&&tuple[2]=='d')
				visit[0][0] = true;
			else
				visit[0][1] = true;
		}
		if(func_num == 2)
		{
			int entityId = (int)x[0];
			double delay = x[1];

			if(entityId<0)
			{
				visit[0][0] = true;
				visit[3][0] = true;
			}
			else
			{
				visit[0][1] = true;
				visit[3][1] = true;
			}
			if(delay<0)
				visit[1][0] = true;
			else
				visit[1][1] = true;
			if(entityId>=999999)
				visit[2][0] = true;
			else
				visit[2][1] = true;
			if(entityId!=1)
				visit[4][0] = true;
			else
				visit[4][1] = true;
		}
		if(func_num == 3)
		{
			int eventTime = (int)x[0];
			int type = (int)x[1];
			int dest = (int)x[2];
			int state = (int)x[3];
			int p = (int)x[4];
			int tag = (int)x[5];
			int src = (int)x[6];

			if(eventTime<50) visit[0][0] = true;
			else visit[0][1] = true;

			if(type==0) visit[1][0] = true;
			else if(type==1) visit[1][1] = true;
			else if(type==2) visit[1][2] = true;
			else visit[1][3] = true;

			if(dest<0) visit[2][0] = true;
			else visit[2][1] = true;

			if(src<0) visit[3][0] = true;
			else visit[3][1] = true;

			if(state==1) visit[4][0] = true;
			else visit[4][1] = true;

			if(p==0||tag==9999) visit[5][0] = true;
			else visit[5][1] = true;
		}
		if(func_num == 4)
		{
			int Direction = (int)x[0];
			String map1 = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			String map2 = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]);

			if(Direction==1)
				visit[0][0] = true;
			else
				visit[0][1] = true;

			if(map1.equals("001")||map1.equals("002")||map1.equals("003"))
				visit[1][0] = true;
			else
				visit[1][1] = true;

			if(map2.equals("004")||map2.equals("005")||map2.equals("006"))
				visit[2][0] = true;
			else
				visit[2][1] = true;
		}
		if(func_num == 5)
		{
			boolean isFinished;
			boolean cloudletCompleted;
			if((int)x[0]==0) isFinished=true;
			else isFinished=false;
			if((int)x[R-1]==0) cloudletCompleted=true;
			else cloudletCompleted=false;
			String cl = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);

			if(isFinished)
				visit[0][0] = true;
			else visit[0][1] = true;

			if(cl.equals("001")||cl.equals("002")||cl.equals("003"))
				visit[1][0] = true;
			else visit[1][1] = true;

			if(cloudletCompleted)
				visit[2][0] = true;
			else visit[2][1] = true;
		}
		if(func_num == 6)
		{
			String edge = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			String pair = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]);
			boolean canSelect;
			if((int)x[6]==0) canSelect = true;
			else canSelect = false;

			if(edge.equals("mod"))
				visit[0][0] = true;
			else visit[0][1] = true;

			if(pair.equals("001")||pair.equals("002")||pair.equals("003"))
				visit[1][0] = true;
			else visit[1][1] = true;

			if(canSelect) visit[2][0] = true;
			else visit[2][1] = true;

			if((int)x[7]==2) visit[3][0] = true;
			else visit[3][1] = true;
		}
		if(func_num == 7)
		{
			String option = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			String extraOption = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]);
			int type = x[6];

			if(option.equals("   "))
				visit[0][0] = true;
			else visit[0][1] = true;

			if(!extraOption.equals("   "))
				visit[1][0] = true;
			else visit[1][1] = true;

			if(extraOption.endsWith(","))
				visit[2][0] = true;
			else visit[2][1] = true;

			if(type == 0) visit[3][0] = true;
			else if(type ==1) visit[3][1] = true;
			else if(type ==2) visit[3][2] = true;
			else if(type ==3) visit[3][3] = true;
			else if(type ==4) visit[3][4] = true;
			else if(type ==5) visit[3][5] = true;
			else if(type ==6) visit[3][6] = true;
			else if(type >6) visit[3][7] = true;
		}
		if(func_num == 8){
			String xmlTags = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			String sentence = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]);

			if(xmlTags.equals("   "))
				visit[0][0] = true;
			else visit[0][1] = true;

			if(sentence.equals("   "))
				visit[1][0] = true;
			else visit[1][1] = true;
		}
		if(func_num == 9){
			Boolean nlSplitting,whitespaceTokenization;
			if(x[0]==0) nlSplitting = true;
			else nlSplitting = false;
			if(x[1]==0) whitespaceTokenization = true;
			else whitespaceTokenization = false;
			String line = String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			String isOneSentence = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6])+String.valueOf((char)x[7]);
			char token,bound1,bound2;
			token = (char)x[8]; bound1 = (char)x[9]; bound2 = (char)x[10];

			if(nlSplitting) visit[0][0] = true;
			else visit[0][1] = true;

			if(whitespaceTokenization) visit[1][0] = true;
			else visit[1][1] = true;

			if(line.equals("/n")) visit[2][0] = true;
			else visit[2][1] = true;

			if(isOneSentence.equals("true")) visit[3][0] = true;
			else visit[3][1] = true;

			if(token==' ') visit[4][0] = true;
			else visit[4][1] = true;

			if(bound1==' ') visit[5][0] = true;
			else visit[5][1] = true;

			if(bound2==' ') visit[6][0] = true;
			else visit[6][1] = true;
		}
		if(func_num == 10){
			String annotation = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			int nThreads = x[3];

			if(annotation.equals("001")||annotation.equals("002")||annotation.equals("003"))
				visit[0][0] = true;
			else
				visit[0][1] = true;

			if(nThreads == 1)
				visit[1][0] = true;
			else
				visit[1][1] = true;
		}
		if(func_num == 11){
			String nerLanguage = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+
					String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]);
			Boolean augment;
			if(x[7]==0) augment = true;
			else augment = false;

			if(nerLanguage.equals("CHINESE"))
				visit[0][0] = true;
			else visit[0][1] = true;

			if(augment) visit[1][0] = true;
			else visit[1][1] = true;
		}
		if(func_num==12){
			String trueCase = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3])+String.valueOf((char)x[4]);
			boolean overwriteText;
			if(x[5]==0) overwriteText = true;
			else overwriteText = false;

			if(trueCase.equals("UPPER")) visit[0][0] = true;
			else visit[0][1] = true;

			if(trueCase.equals("LOWER")) visit[1][0] = true;
			else visit[1][1] = true;

			if(trueCase.equals("INIT_")) visit[2][0] = true;
			else visit[2][1] = true;

			if(trueCase.equals("O    ")) visit[3][0] = true;
			else visit[3][1] = true;

			if(overwriteText) visit[4][0] = true;
			else visit[4][1] = true;
		}
	}

	public static int pathnum(int[] x , int func_num)
	{
		int path = -1;

		if(func_num == 1)
		{
			char tuple[]= new char[R];
			for(int i=0;i<R;i++)
				tuple[i] = (char) x[i];

			if(tuple[0]=='O'&&tuple[1]=='l'&&tuple[2]=='d')
				path = 0;
			else
				path = 1;
		}
		if(func_num == 2)
		{
			int entityId = (int)x[0];
			double delay = x[1];

			if(entityId<0)
			{
				path = 0;
			}else{
				if(delay<0)
				{
					delay = 0;
					if(delay>=999999)
						path = 1;
					else if(entityId<0)
						path = 2;
					else if(entityId!=1)
						path = 3;
					else
						path = 4;

				}
				else{
					if(delay>=999999)
						path = 5;
					else if(entityId<0)
						path = 6;
					else if(entityId!=1)
						path = 7;
					else
						path = 8;
				}
			}
		}
		if(func_num == 3)
		{
			int eventTime = (int) x[0];
			int type = (int)x[1];
			int dest = (int)x[2];
			int state = (int)x[3];
			int p = (int)x[4];
			int tag = (int)x[5];
			int src = (int)x[6];

			if(eventTime<50)
				path = 0;
			else{
				if(type==0)
					path = 1;
				else if(type==3)
					path = 8;
				else if(type==2){
					if(src<0)
						path = 6;
					else
						path = 7;
				}else{
					if(dest<0)
						path = 2;
					else if(state==1){
						path = 5;
					}else{
						if(p==0||tag==9999)
							path = 3;
						else
							path = 4;
					}
				}
			}
		}
		if(func_num == 4)
		{
			int Direction = (int)x[0];
			String map1 = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			String map2 = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]);

			if(Direction==1){
				if(map1.equals("001")||map1.equals("002")||map1.equals("003")){
					if(map2.equals("004")||map2.equals("005")||map2.equals("006"))
						path = 0;
					else
						path = 1;
				}else{
					if(map2.equals("004")||map2.equals("005")||map2.equals("006"))
						path = 2;
					else
						path = 3;
				}
			}
			else path = 4;
		}
		if(func_num == 5)
		{
			boolean isFinished;
			boolean cloudletCompleted;
			if((int)x[0]==0) isFinished=true;
			else isFinished=false;
			if((int)x[R-1]==0) cloudletCompleted=true;
			else cloudletCompleted=false;
			String cl = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);

			if(isFinished)
			{
				if(cl.equals("001")||cl.equals("002")||cl.equals("003"))
					if(cloudletCompleted)
						path = 0;
					else
						path = 1;
				else
				if(cloudletCompleted)
					path = 2;
				else
					path = 3;
			}else
			if(cloudletCompleted)
				path = 4;
			else
				path = 5;
		}
		if(func_num == 6)
		{
			String edge = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			String pair = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]);
			boolean canSelect;
			if((int)x[6]==0) canSelect = true;
			else canSelect = false;

			if(edge.equals("mod"))
			{
				if(pair.equals("001")||pair.equals("002")||pair.equals("003")){
					if(canSelect){
						if((int)x[7]==2)
							path = 0;
						else
							path = 1;
					}else
						path = 2;
				}else{
					if(canSelect){
						if((int)x[7]==2)
							path = 3;
						else
							path = 4;
					}else
						path = 5;
				}
			}else
				path = 6;
		}

		if(func_num == 7)
		{
			String option = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			String extraOption = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]);
			int type = x[6];

			if(option.equals("   "))
			{
				if(!extraOption.equals("   ")){
					if(extraOption.endsWith(","))
					{//"0000"-"0007"
						if(type == 0) path = 0;
						else if(type == 1) path = 1;
						else if(type == 2) path = 2;
						else if(type == 3) path = 3;
						else if(type == 4) path = 4;
						else if(type == 5) path = 5;
						else if(type == 6) path = 6;
						else if(type >6) path = 7;
					}
					else{//"0010"-"0017"
						if(type == 0) path = 8;
						else if(type == 1) path = 9;
						else if(type == 2) path = 10;
						else if(type == 3) path = 11;
						else if(type == 4) path = 12;
						else if(type == 5) path = 13;
						else if(type == 6) path = 14;
						else if(type >6) path = 15;
					}
				}else{//"01 0"-"01 7"
					if(type == 0) path = 16;
					else if(type == 1) path = 17;
					else if(type == 2) path = 18;
					else if(type == 3) path = 19;
					else if(type == 4) path = 20;
					else if(type == 5) path = 21;
					else if(type == 6) path = 22;
					else if(type >6) path = 23;
				}
			}else{
				if(!extraOption.equals("   ")){
					if(extraOption.endsWith(","))
					{//"1000"-"1007"
						if(type == 0) path = 24;
						else if(type == 1) path = 25;
						else if(type == 2) path = 26;
						else if(type == 3) path = 27;
						else if(type == 4) path = 28;
						else if(type == 5) path = 29;
						else if(type == 6) path = 30;
						else if(type >6) path = 31;
					}
					else{//"1010"-"1017"
						if(type == 0) path = 32;
						else if(type == 1) path = 33;
						else if(type == 2) path = 34;
						else if(type == 3) path = 35;
						else if(type == 4) path = 36;
						else if(type == 5) path = 37;
						else if(type == 6) path = 38;
						else if(type >6) path = 39;
					}
				}else{//"11 0"-"11 7"
					if(type == 0) path = 40;
					else if(type == 1) path = 41;
					else if(type == 2) path = 42;
					else if(type == 3) path = 43;
					else if(type == 4) path = 44;
					else if(type == 5) path = 45;
					else if(type == 6) path = 46;
					else if(type >6) path = 47;
				}
			}

		}
		if(func_num == 8){
			String xmlTags = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			String sentence = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]);

			if(xmlTags.equals("   "))
			{
				if(sentence.equals("   "))
					path = 0;
				else
					path = 1;
			}else
				path = 2;
		}
		if(func_num == 9){
			Boolean nlSplitting,whitespaceTokenization;
			if(x[0]==0) nlSplitting = true;
			else nlSplitting = false;
			if(x[1]==0) whitespaceTokenization = true;
			else whitespaceTokenization = false;
			String line = String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			String isOneSentence = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6])+String.valueOf((char)x[7]);
			char token,bound1,bound2;
			token = (char)x[8]; bound1 = (char)x[9]; bound2 = (char)x[10];

			if(nlSplitting){
				if(whitespaceTokenization){
					if((line.equals("/n")))
						path = 0;
					else
						path = 1;
				}
				else
					path = 2;
			}
			else{
				if(isOneSentence.equals("true"))
					path = 3;
				else{
					if(token==' '){
						if(bound1==' '){
							if(bound2==' ')
								path = 4;
							else
								path = 5;
						}else{
							if(bound2==' ')
								path = 6;
							else
								path = 7;
						}
					}else{
						if(bound1==' '){
							if(bound2==' ')
								path = 8;
							else
								path = 9;
						}else{
							if(bound2==' ')
								path = 10;
							else
								path = 11;
						}
					}
				}
			}
		}
		if(func_num == 10)
		{
			String annotation = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			int nThreads = x[3];

			if(annotation.equals("001")||annotation.equals("002")||annotation.equals("003"))
			{
				if(nThreads == 1)
					path =0;
				else
					path = 1;
			}else
				path =2;
		}
		if(func_num == 11){
			String nerLanguage = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+
					String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]);
			Boolean augment;
			if(x[7]==0) augment = true;
			else augment = false;

			if(nerLanguage.equals("CHINESE")){
				if(augment)
					path = 0;
				else
					path = 1;
			}else{
				if(augment)
					path = 2;
				else
					path = 3;
			}
		}
		if(func_num==12){
			String trueCase = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3])+String.valueOf((char)x[4]);
			boolean overwriteText;
			if(x[5]==0) overwriteText = true;
			else overwriteText = false;

			if(trueCase.equals("UPPER")){
				if(overwriteText)
					path = 0;
				else
					path = 1;
			}else if(trueCase.equals("LOWER")){
				if(overwriteText)
					path = 2;
				else
					path = 3;
			}else if(trueCase.equals("INIT_")){
				if(overwriteText)
					path = 4;
				else
					path = 5;
			}else if(trueCase.equals("O    ")){
				if(overwriteText)
					path = 6;
				else
					path = 7;
			}else{
				if(overwriteText)
					path = 8;
				else
					path = 9;
			}
		}

		return path;
	}

	public static double benchmarkfunction (int[] x , int func_num, int path_num)
	{
		double[] fit = new double[NODENUM] ;
		double[][] f = new double[NODENUM][BRANCH];
		double Fitness = 0 ;

		if(func_num == 1)
		{
			char tuple[]= new char[R];
			for(int i=0;i<R;i++)
				tuple[i] = (char) x[i];

			double v1=0;

			if(tuple[0]=='O'&&tuple[1]=='l'&&tuple[2]=='d')
				v1 = 0;
			else
				v1 = Math.abs(tuple[0]-'O')+K+ Math.abs(tuple[1]-'l')+K + Math.abs(tuple[2]-'d')+K;
			f[0][0] = v1;

			if(tuple[0]!='O'||tuple[1]!='l'||tuple[2]!='d')
				v1 = 0;
			else{
				v1 = Math.min(Math.abs(tuple[0]-'O'), Math.abs(tuple[1]-'l'));
				v1 = Math.min(v1, Math.abs(tuple[2]-'d'));
			}
			f[0][1] = v1;
		}
		if(func_num == 2)
		{
			int entityId = (int)x[0];
			double delay = x[1];

			double v1=0,v2=0;

			//分支节点1
			if(entityId<0) v1=0;
			else v1 = entityId+K;
			f[0][0] = v1;

			if(entityId>=0) v2=0;
			else v2 = 0-entityId+K;
			f[0][1] = v2;

			//分支节点2
			if(delay<0) {v1=0;delay=0;}
			else v1 = delay+K;
			f[1][0] = v1;

			if(delay>=0) v2=0;
			else v2 = -delay+K;
			f[1][1] = v2;

			//分支节点3
			if(delay>=999999) v1=0;
			else v1 = 999999-delay+K;
			f[2][0] = v1;

			if(delay<999999) v2=0;
			else v2 = delay-999999+K;
			f[2][1] = v2;

			//分支节点4
			if(entityId<0) v1=0;
			else v1 = entityId+K;
			f[3][0] = v1;

			if(entityId>=0) v2 = 0;
			else v2 = 0-entityId+K;
			f[3][1] = v2;

			//分支节点5
			if(entityId!=1) v1=0;
			else v1 = K;
			f[4][0] = v1;

			if(entityId==1) v2 = 0;
			else v2 = Math.abs(entityId+1)+K;
			f[4][1] = v2;
		}
		if(func_num == 3)
		{
			int eventTime = (int)x[0];
			int type = (int)x[1];
			int dest = (int)x[2];
			int state = (int)x[3];
			int p = (int)x[4];
			int tag = (int)x[5];
			int src = (int)x[6];

			double v1=0;
			//分支节点0
			if(eventTime<50) v1=0;
			else v1 = eventTime-50+K;
			f[0][0] = v1;

			if(eventTime>=50) v1=0;
			else v1 = 50-eventTime+K;
			f[0][1] = v1;

			//分支节点1
			if(type==0)	f[1][0] = 0;
			else	f[1][0] = Math.abs(type-0)+K;
			if(type==1)	f[1][1] = 0;
			else	f[1][1] = Math.abs(type-1)+K;
			if(type==2)	f[1][2] = 0;
			else	f[1][2] = Math.abs(type-2)+K;
			if(type==3)	f[1][3] = 0;
			else	f[1][3] = Math.abs(type-3)+K;


			//分支节点2
			if(dest<0) v1=0;
			else v1 = dest+K;
			f[2][0] = v1;

			if(dest>=0) v1=0;
			else v1=-dest+K;
			f[2][1] = v1;

			//分支节点3
			if(src<0) v1=0;
			else v1=src+K;
			f[3][0] = v1;

			if(src>=0) v1=0;
			else v1= -src+K;
			f[3][1] = v1;

			//分支节点4
			if(state==1)v1=0;
			else v1=Math.abs(state-1)+K;
			f[4][0] = v1;

			if(state==0) v1=0;
			else v1= Math.abs(state)+K;
			f[4][1] = v1;

			//分支节点5
			if(p==0||tag==9999)v1=0;
			else v1=Math.min(Math.abs(p-0)+K,Math.abs(tag-9999)+K);
			f[5][0] = v1;

			if(p!=0&&tag!=9999)v1=0;
			else v1= 2*K;
			f[5][1] = v1;
		}
		if(func_num == 4)
		{
			int Direction = (int)x[0];
			String map1 = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			String map2 = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]);

			double v1=0;

			//分支节点1
			if(Direction==1) v1=0;
			else v1 = Math.abs(Direction-1)+K;
			f[0][0] = v1;

			if(Direction!=1) v1=0;
			else v1 = K;
			f[0][1] = v1;

			//分支节点2
			if(map1.equals("001")||map1.equals("002")||map1.equals("003"))
				v1=0;
			else
//				v1 = Math.abs((char)x[1]-'0') + Math.abs((char)x[2]-'0') + Math.min(Math.min((char)x[3]-'1', (char)x[3]-'2'), (char)x[3]-'3');
			{
				v1 = Math.min(Math.abs((char)x[3]-'1')+K, Math.abs((char)x[3]-'2')+K);
				v1 = Math.min(v1, Math.abs((char)x[3]-'3')+K);
				v1 = Math.abs((char)x[1]-'0') + Math.abs((char)x[2]-'0') + v1;
			}
			f[1][0] = v1;

			if(!map1.equals("001")&&!map1.equals("002")&&!map1.equals("003"))
				v1 = 0;
			else
				v1 = 3*K;
			f[1][1] = v1;

			//分支节点3
			if(map2.equals("004")||map2.equals("005")||map2.equals("006"))
				v1=0;
			else
//				v1 = Math.abs((char)x[4]-'0') + Math.abs((char)x[5]-'0') + Math.min(Math.min((char)x[6]-'4', (char)x[6]-'5'), (char)x[6]-'6');
			{
				v1 = Math.min(Math.abs((char)x[6]-'4')+K, Math.abs((char)x[6]-'5')+K);
				v1 = Math.min(v1, Math.abs((char)x[6]-'6')+K);
				v1 = Math.abs((char)x[4]-'0') + Math.abs((char)x[5]-'0') + v1;
			}
			f[2][0] = v1;

			if(!map2.equals("004")&&!map2.equals("005")&&!map2.equals("006"))
				v1 = 0;
			else
				v1 = 3*K;
			f[2][1] = v1;
		}
		if(func_num == 5)
		{
			String cl = String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			int p1 = (int)x[0], p2 = (int)x[R-1];

			double v1=0;
			//分支节点1
			if(p1==0) v1=0;
			else v1=Math.abs(p1)+K;
			f[0][0] = v1;

			if(p1==1) v1=0;
			else v1 = Math.abs(p1-1)+K;
			f[0][1] = v1;

			//分支节点2
			if(cl.equals("001")||cl.equals("002")||cl.equals("003"))
				v1=0;
			else
//				v1=Math.abs((char)x[1]-'0') + Math.abs((char)x[2]-'0') + Math.min(Math.min((char)x[3]-'1', (char)x[3]-'2'), (char)x[3]-'3');
			{
				v1 = Math.min(Math.abs((char)x[3]-'1')+K, Math.abs((char)x[3]-'2')+K);
				v1 = Math.min(v1, Math.abs((char)x[3]-'3')+K);
				v1 = Math.abs((char)x[1]-'0') + Math.abs((char)x[2]-'0') + v1;
			}
			f[1][0] = v1;

			if(!cl.equals("001")&&!cl.equals("002")&&!cl.equals("003"))
				v1 = 0;
			else
				v1 = 3*K;
			f[1][1] = v1;

			//分支节点3
			if(p2==0) v1=0;
			else v1=Math.abs(p2)+K;
			f[2][0] = v1;

			if(p2==1) v1=0;
			else v1 = Math.abs(p2-1)+K;
			f[2][1] = v1;
		}
		if(func_num == 6)
		{
			String edge = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			String pair = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]);

			double v1=0;
			//分支节点1
			if(edge.equals("mod")) v1=0;
			else v1 = Math.abs((char)x[0]-'m')+Math.abs((char)x[1]-'o')+Math.abs((char)x[2]-'d')+3*K;
			f[0][0] = v1;

			if(!edge.equals("mod")) v1=0;
			else v1 = K;
			f[0][1] = v1;

			//分支节点2
			if(pair.equals("001")||pair.equals("002")||pair.equals("003"))
				v1 = 0;
			else
//				v1 = Math.abs((char)x[3]-'0') + Math.abs((char)x[4]-'0') + Math.min(Math.min((char)x[5]-'1', (char)x[5]-'2'), (char)x[5]-'3');
			{
				v1 = Math.min(Math.abs((char)x[5]-'1')+K, Math.abs((char)x[5]-'2')+K);
				v1 = Math.min(v1, Math.abs((char)x[5]-'3')+K);
				v1 = Math.abs((char)x[3]-'0') + Math.abs((char)x[4]-'0') + v1;
			}
			f[1][0] = v1;

			if(!pair.equals("001")&&!pair.equals("002")||!pair.equals("003"))
				v1 = 0;
			else
				v1 = 3*K;
			f[1][1] = v1;
			//分支节点3
			if((int)x[6]==0) v1 = 0;
			else v1 = Math.abs((int)x[6])+K;
			f[2][0] = v1;

			if((int)x[6]==1) v1 =0;
			else v1 = Math.abs((int)x[6]-1)+K;
			f[2][1] = v1;
			//分支节点4
			if((int)x[7]==2) v1=0;
			else v1 = Math.abs((int)x[7]-2)+K;
			f[3][0] = v1;

			if((int)x[7]!=2) v1=0;
			else v1 = K;
			f[3][1] = K;
		}


		if(func_num == 7)
		{
			String option = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			String extraOption = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]);
			int type = x[6];

			double v1=0;
			//分支节点1
			if(option.equals("   ")) v1 = 0;
			else v1 = Math.abs((char)x[0]-' ')+Math.abs((char)x[1]-' ')+Math.abs((char)x[2]-' ')+3*K;
			f[0][0] = v1;

			if(!option.equals("   ")) v1 =0;
			else v1 = K;
			f[0][1] = v1;

			//分支节点2
			if(!extraOption.equals("   ")) v1=0;
			else v1 = K;
			f[1][0] = v1;

			if(extraOption.equals("   ")) v1=0;
			else v1 = Math.abs((char)x[3]-' ')+Math.abs((char)x[4]-' ')+Math.abs((char)x[5]-' ')+3*K;
			f[1][1] = v1;

			//分支节点3
			if(x[5]==',') v1=0;
			else v1 = Math.abs((char)x[5] - ',') + K;
			f[2][0] = v1;

			if(x[5]!=',') v1 = 0;
			else v1 = K;
			f[2][1] = v1;

			//分支节点4
			if(type==0)	f[3][0] = 0;
			else	f[3][0] = Math.abs(type-0)+K;
			if(type==1)	f[3][1] = 0;
			else	f[3][1] = Math.abs(type-1)+K;
			if(type==2)	f[3][2] = 0;
			else	f[3][2] = Math.abs(type-2)+K;
			if(type==3)	f[3][3] = 0;
			else	f[3][3] = Math.abs(type-3)+K;
			if(type==4)	f[3][4] = 0;
			else	f[3][4] = Math.abs(type-4)+K;
			if(type==5)	f[3][5] = 0;
			else	f[3][5] = Math.abs(type-5)+K;
			if(type==6)	f[3][6] = 0;
			else	f[3][6] = Math.abs(type-6)+K;
		}
		if(func_num == 8){
			String xmlTags = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			String sentence = String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5]);

			double v1=0;
			//分支节点1
			if(xmlTags.equals("   ")) v1 = 0;
			else v1 = Math.abs((char)x[0]-' ')+Math.abs((char)x[1]-' ')+Math.abs((char)x[2]-' ')+3*K;
			f[0][0] = v1;

			if(!xmlTags.equals("   ")) v1 =0;
			else v1 = K;
			f[0][1] = v1;

			//分支节点2
			if(sentence.equals("   ")) v1 = 0;
			else v1 = Math.abs((char)x[3]-' ')+Math.abs((char)x[4]-' ')+Math.abs((char)x[5]-' ')+3*K;
			f[1][0] = v1;

			if(!sentence.equals("   ")) v1 =0;
			else v1 = K;
			f[1][1] = v1;
		}
		if(func_num == 9){
			boolean nlSplitting,whitespaceTokenization;
			if(x[0]==0) nlSplitting = true;
			else nlSplitting = false;
			if(x[1]==0) whitespaceTokenization = true;
			else whitespaceTokenization = false;
			String line = String.valueOf((char)x[2])+String.valueOf((char)x[3]);
			String isOneSentence = String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6])+String.valueOf((char)x[7]);
			char token,bound1,bound2;
			token = (char)x[8]; bound1 = (char)x[9]; bound2 = (char)x[10];

			double v1=0;
			//分支节点1
			if(nlSplitting) v1 = 0;
			else v1 = K;
			f[0][0] = v1;

			if(!nlSplitting) v1 = 0;
			else v1 = K;
			f[0][1] = v1;

			//分支节点2
			if(whitespaceTokenization) v1 = 0;
			else v1 = K;
			f[1][0] = v1;

			if(!whitespaceTokenization) v1 = 0;
			else v1 = K;
			f[1][1] = v1;

			//分支节点3
			if(line.equals("/n")) v1 = 0;
			else v1 = Math.abs((char)x[2]-'/')+Math.abs((char)x[3]-'n')+2*K;
			f[2][0] = v1;

			if(!line.equals("\n")) v1 = 0;
			else v1 = K;
			f[2][1] = v1;

			//分支节点4
			if(isOneSentence.equals("true")) v1 = 0;
			else v1 = Math.abs((char)x[4]-'t')+Math.abs((char)x[5]-'r')+Math.abs((char)x[6]-'u')+Math.abs((char)x[7]-'e')+4*K;
			f[3][0] = v1;

			if(!isOneSentence.equals("true")) v1 = 0;
			else v1 = K;
			f[3][1] = v1;

			//分支节点5
			if(token==' ') v1 = 0;
			else v1 = Math.abs(token-' ')+K;
			f[4][0] = v1;

			if(token!=' ') v1 = 0;
			else v1 = K;
			f[4][1] = v1;

			//分支节点6
			if(bound1==' ') v1 = 0;
			else v1 = Math.abs(bound1-' ')+K;
			f[5][0] = v1;

			if(bound1!=' ') v1 = 0;
			else v1 = K;
			f[5][1] = v1;

			//分支节点7
			if(bound2==' ') v1 = 0;
			else v1 = Math.abs(bound2-' ')+K;
			f[6][0] = v1;

			if(bound2!=' ') v1 = 0;
			else v1 = K;
			f[6][1] = v1;
		}
		if(func_num == 10){
			String annotation = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2]);
			int nThreads = x[3];

			double v1;
			//分支节点1
			if(annotation.equals("001")||annotation.equals("002")||annotation.equals("003")) v1 = 0;
			else {
				v1 = Math.min(Math.abs((char)x[2]-'1')+K, Math.abs((char)x[2]-'2')+K);
				v1 = Math.min(v1, Math.abs((char)x[2]-'3')+K);
				v1 = Math.abs((char)x[0]-'0') + Math.abs((char)x[1]-'0') + v1;
			}
			f[0][0] = v1;

			if(!annotation.equals("001")&&!annotation.equals("002")||!annotation.equals("003"))
				v1 = 0;
			else
				v1 = 3*K;
			f[0][1] = v1;
			//分支节点2
			if(nThreads == 1) v1 = 0;
			else v1 = Math.abs(nThreads-1)+K;
			f[1][0] = v1;

			if(nThreads != 1) v1=0;
			else v1 = K;
			f[1][1] = v1;
		}
		if(func_num == 11){
			String nerLanguage = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+
					String.valueOf((char)x[3])+String.valueOf((char)x[4])+String.valueOf((char)x[5])+String.valueOf((char)x[6]);
			boolean augment;
			if(x[7]==0) augment = true;
			else augment = false;

			double v1;
			//分支节点1
			if(nerLanguage.equals("CHINESE")) v1 = 0;
			else v1 = Math.abs((char)x[0]-'C')+Math.abs((char)x[1]-'H')+Math.abs((char)x[2]-'I')+Math.abs((char)x[3]-'N')+
					Math.abs((char)x[4]-'E')+Math.abs((char)x[5]-'S')+Math.abs((char)x[6]-'E')+7*K;
			f[0][0] = v1;

			if(!nerLanguage.equals("CHINESE")) v1 = 0;
			else v1 = 7*K;
			f[0][1] = v1;

			//分支节点2
			if(augment) v1=0;
			else v1 = K;
			f[1][0] = v1;

			if(!augment) v1=0;
			else v1 = K;
			f[1][1] = v1;
		}
		if(func_num == 12){
			String trueCase = String.valueOf((char)x[0])+String.valueOf((char)x[1])+String.valueOf((char)x[2])+String.valueOf((char)x[3])+String.valueOf((char)x[4]);
			boolean overwriteText;
			if(x[5]==0) overwriteText = true;
			else overwriteText = false;

//			double v1;
			//分支节点1
			if(trueCase.equals("UPPER"))
				f[0][0] = 0;
			else
				f[0][0] = Math.abs((char)x[0]-'U')+Math.abs((char)x[1]-'P')+Math.abs((char)x[2]-'P')+Math.abs((char)x[3]-'E')+Math.abs((char)x[4]-'R')+5*K;

			if(!trueCase.equals("UPPER"))
				f[0][1] = 0;
			else
				f[0][1] = 5*K;

			//分支节点2
			if(trueCase.equals("LOWER"))
				f[1][0] = 0;
			else
				f[1][0] = Math.abs((char)x[0]-'L')+Math.abs((char)x[1]-'O')+Math.abs((char)x[2]-'W')+Math.abs((char)x[3]-'E')+Math.abs((char)x[4]-'R')+5*K;

			if(!trueCase.equals("LOWER"))
				f[1][1] = 0;
			else
				f[1][1] = 5*K;

			//分支节点3
			if(trueCase.equals("INIT_"))
				f[2][0] = 0;
			else
				f[2][0] = Math.abs((char)x[0]-'I')+Math.abs((char)x[1]-'N')+Math.abs((char)x[2]-'I')+Math.abs((char)x[3]-'T')+Math.abs((char)x[4]-'_')+5*K;

			if(!trueCase.equals("INIT_"))
				f[2][1] = 0;
			else
				f[2][1] = 5*K;

			//分支节点4
			if(trueCase.equals("O    "))
				f[3][0] = 0;
			else
				f[3][0] = Math.abs((char)x[0]-'O')+Math.abs((char)x[1]-' ')+Math.abs((char)x[2]-' ')+Math.abs((char)x[3]-' ')+Math.abs((char)x[4]-' ')+5*K;

			if(!trueCase.equals("O   "))
				f[3][1] = 0;
			else
				f[3][1] = 5*K;

			//分支节点5
			if(overwriteText) {f[4][0] = 0; f[4][1] = K;}
			else {f[4][0] = K; f[4][1] = 0;}
		}
		if(path_num == -1)          //没有目标路径的情况
		{
			for(int k = 0 ; k < NODENUM ; k++)
			{
				if(visit[k][0] && visit[k][1])
					fit[k] = 0 ;
				else if(visit[k][0] && (!visit[k][1]))
					fit[k] = 1/(f[k][1] + alpha) ;
				else if((!visit[k][0]) && visit[k][1])
					fit[k] = 1/(f[k][0] + alpha) ;
				else
					fit[k] = 1/alpha ;
			}
		}
		else {// 存在目标路径的情况
			for (int k = 0; k < NODENUM; k++) {
				if (PATH[path_num].charAt(k) == ' ')
					fit[k] = 0;
				else
					fit[k] = 1 / (f[k][PATH[path_num].charAt(k) - '0'] + alpha);

			}

		}
		for(int i=0;i<NODENUM;i++)
			Fitness += fit[i];

		return Fitness;
	}


	public static int[] selectsort (double[] fitness)
	{
		int[] sortnum = new int[fitness.length] ;
		for(int i =0 ; i < fitness.length ; i++)
			sortnum[i]=i;
		int max = 0;
		double tmp = 0;
		int tmp2 = 0 ;
		for(int i=0;i<fitness.length;i++)
		{
			max = i;

			for(int j=i+1;j<fitness.length;j++)
			{
				if(fitness[max]<fitness[j])
					max = j;
			}
			if(i!=max)
			{
				tmp = fitness[i];
				fitness[i] = fitness[max];
				fitness[max] = tmp;

				tmp2 = sortnum[i];
				sortnum[i] = sortnum[max];
				sortnum[max] = tmp2;
			}
		}
		return sortnum ;
	}


	static double getAverage(int[] array , int num){
		int sum = 0;
		for(int i = 0;i < num;i++){
			sum += array[i];
		}
		return (double)(sum / num);
	}


	static double getStandardDevition(int[] array , int num){
		double sum = 0;
		for(int i = 0;i < num;i++){
			sum += Math.sqrt(((double)array[i] -getAverage(array, num)) * (array[i] -getAverage(array, num)));
		}
		return (sum / (num - 1));
	}

	static int getBestIndex(double[]array)
	{
		int index;
		index = 0;
		for(int i=1;i<array.length;i++)
			if(array[i]>array[index])
				index = i;
		return index;
	}
	static int[] getIndex(int lb, int ub, int best, double step) {
		int[] index = new int[step_length+1];
		int temp = -1;
		if (best + step_length/2 * step <= ub && best - step_length/2 * step >= lb) {
			for (int i = 0; i < step_length+1; i++)
				index[i] = (int) (step * (i - step_length/2) + best);
		} else {
			for (int i = 1; i < step_length+1; i++)
				if (best - i * step < lb) {
					temp = i;
					break;
				}
			if (temp == -1) {
				for (int i = 0; i < step_length+1; i++)
					index[i] = (int) (best - step * i) ;
			} else {
				for (int i = 0; i < step_length+1; i++)
					index[i] = (int) (best + step * (-temp + i + 1));
			}
		}

		return index;
	}

}
