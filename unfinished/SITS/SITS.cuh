#ifndef SITS_CUH
#define SITS_CUH
#include "../control.cuh"
#include "../common.cuh"
enum SITS_MODE
{
    ONLY_SELECT_MODE_FORCE = 0,
    SIMPLE_SITS_MODE = 1,
    ONLY_SELECT_MODE_FORCE_ENERGY = 2,
	CLASSICAL_SITS_MODE = 3,

};
struct SITS
{
private:
	FILE *sits_ene_record_out = NULL;
	struct FC_BALL_INFORMATION
	{
		float move_length = 0.01;//simpleSITS中fcball随机游走的最大步长
		float fc_max = 1.2;//游走的上限，对应最低的温度T正比于1/fc_ball
		float fc_min = 0.5;//游走的下限，对应最高的温度

		int random_seed = 0;//随机游走的初始种子，可能和其他程序的种子冲突
		float *fc_pdf = NULL;//离散的fc概率密度记录列表，用于控制fc的概率分布，cpu上存储
		int grid_numbers = 1000;//离散列表的格子数目

		float current_fc_probability=0.01;//初始的概率密度,要非0
		float get_fc_probability(float pos);//获得此时fc_ball的值pos对应的概率密度
		

		int is_constant_fc_ball=0;//记录是否是固定fcball值进行模拟（用于体系初测）
		float constant_fc_ball = 1.0;//固定的fcball值
	}simple_info;
    void fc_ball_random_walk();//simple mode里根据上面几个参数进行fc_ball的一次随机移动
	void SITS_Classical_Update_Info(int steps);//classical info中需要迭代更新Nk
public:
	struct CLASSICAL_SITS_INFORMATION
	{
	public:
		//暂时变量
		int record_count = 0;     //记录次数
		int reset = 1;  //record的时候，第一次和后面公式不一样，这个变量是拿来控制这个的
	
		//控制变量
		int record_interval = 1;  //每隔1步记录一次能量
		int update_interval = 100; //每隔100步更新一次nk
		int constant_nk = 0;         //sits是否迭代更新nk
		int k_numbers;           //划分多少个格子
		float beta0;             //本身温度对应的beta
		//文件
		FILE *nk_traj_file; //记录nk变化的文件
		char nk_rest_file[256]; //记录最后一帧nk的文件
		FILE *norm_traj_file; //记录log_norm变化的文件
		char norm_rest_file[256]; //记录最后一帧log_norm的文件
		float *log_nk_recorded_cpu;  //在cpu端的记录值
		float *log_norm_recorded_cpu; //在cpuu端的记录值

		//计算时，可以对fc_ball直接修正，+ fb_shift进行调节，
		float fb_shift;
		//也可以对进行修正，使加强计算能量时值为 energy_multiple * 原始能量 + energy_shift;
		float energy_multiple;
		float energy_shift;

		////原理和物理量见下面两篇文献
		////A selective integrated tempering method
		////Self-adaptive enhanced sampling in the energy and trajectory spaces : Accelerated thermodynamics and kinetic calculations

		float *beta_k;           
		float *NkExpBetakU;      
		float *Nk;            
		float *sum_a;
		float *sum_b;
		float *d_fc_ball;
		//xyj的cpp变量名-ylj的F90变量名-文献对应
		//ene_recorded - vshift - ene
		//gf - gf - log( n_k * exp(-beta_k * ene) )
		//gfsum - gfsum - log( Sum_(k=1)^N ( log( n_k * exp(-beta_k * ene) ) ) )
		//log_weight - rb - log of the weighting function
		//log_mk_inverse - ratio - log(m_k^-1)
		//log_norm_old - normlold - W(j-1)
		//log_norm - norml - W(j)
		//log_pk - rbfb - log(p_k)
		//log_nk_inverse - pratio - log(n_k^-1)
		//log_nk - fb - log(n_k)
		float *ene_recorded;
		float *gf;
		float *gfsum;
		float *log_weight;
		float *log_mk_inverse;
		float *log_norm_old;
		float *log_norm;
		float *log_pk;
		float *log_nk_inverse;
		float *log_nk;

		void Export_Restart_Information_To_File();
	}classical_info;
public:
	struct SITS_INFORMATION
	{
		int sits_mode = 0;//选择SITS模式
		int atom_numbers = 0;//体系的所有原子数，但不一定要和MD_INFORMATION一样，根据需要
		int protein_atom_numbers = 0;//进行SITS增强的原子数，原子序数一定要是[0-protein_atom_numbers)里的，因此需要预先排序
		int max_neighbor_numbers = 800;//这个数一般与neighbor_list里面的一样
		float pwwp_enhance_factor = 0.5; //相互作用能部分有多少加强
        float fc_ball = 1.;//SITS的强化因子
	}info;

	
	float *d_total_pp_atom_energy = NULL;//AA总能量
	float *d_total_ww_atom_energy = NULL;//BB总能量
	float h_total_pp_atom_energy = 0;//AA总能量
	float h_total_ww_atom_energy = 0;//BB总能量    
    
	float *d_protein_water_atom_energy = NULL;//用于记录每个原子的AB分能量
	float *d_total_protein_water_atom_energy = NULL;//AB总能量
    float h_total_protein_water_atom_energy = 0;//AB总能量
	VECTOR *protein_water_frc = NULL;//用于记录AB两类原子交叉项作用力
	ATOM_GROUP *d_nl_ppww = NULL;//A的A邻居与B的B邻居
	ATOM_GROUP *d_nl_pwwp = NULL;//A的B邻居与B的A邻居

	//初始化
	void Initial_SITS(CONTROLLER *controller, int atom_numbers);
	//清除上一步算的力、能量，明确下一步是否需要计算能量
	void Prepare_For_Calculate_Force(int *need_atom_energy, int isPrintStep);
	
	//根据统计热力学原理，增强的力只需要和能量匹配就好，因此不需要对所有的相互作用进行增强
	//20201217，目前主要增强的力是bond,angle,dihedral,LJ，PME_Direct,nb_14。其他部分不增强（不影响结果正确性）
	//在计算各种SITS相关的能量之前先清空一次
	void Clear_SITS_Energy();
	


	//在上下两个函数中间就可以填入所有需要增强的能量计算子模块的原子能量计算函数
	//在计算完各SITS分能量后进行能量求和
	void Calculate_Total_SITS_Energy(float *d_atom_energy);
    
	//改变frc，使得需要增强的frc被选择性增强，由于共用的frc，因此这步frc增强需要放到刚好计算完所有要增强的frc的函数下方，而避免增强不应该增强的子模块frc
	//因此，体系的所有可能的frc需要先算待增强的frc，插入此函数，再算不应增强的frc
	void SITS_Enhanced_Force(int steps, VECTOR *frc);

	//将信息打印
	void Print();

	//程序结束后，释放SITS内的各个指针
	void Clear_SITS();
};
#endif //SITS_CUH(sits.cuh)
