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
		float move_length = 0.01;//simpleSITS��fcball������ߵ���󲽳�
		float fc_max = 1.2;//���ߵ����ޣ���Ӧ��͵��¶�T������1/fc_ball
		float fc_min = 0.5;//���ߵ����ޣ���Ӧ��ߵ��¶�

		int random_seed = 0;//������ߵĳ�ʼ���ӣ����ܺ�������������ӳ�ͻ
		float *fc_pdf = NULL;//��ɢ��fc�����ܶȼ�¼�б����ڿ���fc�ĸ��ʷֲ���cpu�ϴ洢
		int grid_numbers = 1000;//��ɢ�б�ĸ�����Ŀ

		float current_fc_probability=0.01;//��ʼ�ĸ����ܶ�,Ҫ��0
		float get_fc_probability(float pos);//��ô�ʱfc_ball��ֵpos��Ӧ�ĸ����ܶ�
		

		int is_constant_fc_ball=0;//��¼�Ƿ��ǹ̶�fcballֵ����ģ�⣨������ϵ���⣩
		float constant_fc_ball = 1.0;//�̶���fcballֵ
	}simple_info;
    void fc_ball_random_walk();//simple mode��������漸����������fc_ball��һ������ƶ�
	void SITS_Classical_Update_Info(int steps);//classical info����Ҫ��������Nk
public:
	struct CLASSICAL_SITS_INFORMATION
	{
	public:
		//��ʱ����
		int record_count = 0;     //��¼����
		int reset = 1;  //record��ʱ�򣬵�һ�κͺ��湫ʽ��һ��������������������������
	
		//���Ʊ���
		int record_interval = 1;  //ÿ��1����¼һ������
		int update_interval = 100; //ÿ��100������һ��nk
		int constant_nk = 0;         //sits�Ƿ��������nk
		int k_numbers;           //���ֶ��ٸ�����
		float beta0;             //�����¶ȶ�Ӧ��beta
		//�ļ�
		FILE *nk_traj_file; //��¼nk�仯���ļ�
		char nk_rest_file[256]; //��¼���һ֡nk���ļ�
		FILE *norm_traj_file; //��¼log_norm�仯���ļ�
		char norm_rest_file[256]; //��¼���һ֡log_norm���ļ�
		float *log_nk_recorded_cpu;  //��cpu�˵ļ�¼ֵ
		float *log_norm_recorded_cpu; //��cpuu�˵ļ�¼ֵ

		//����ʱ�����Զ�fc_ballֱ��������+ fb_shift���е��ڣ�
		float fb_shift;
		//Ҳ���ԶԽ���������ʹ��ǿ��������ʱֵΪ energy_multiple * ԭʼ���� + energy_shift;
		float energy_multiple;
		float energy_shift;

		////ԭ�����������������ƪ����
		////A selective integrated tempering method
		////Self-adaptive enhanced sampling in the energy and trajectory spaces : Accelerated thermodynamics and kinetic calculations

		float *beta_k;           
		float *NkExpBetakU;      
		float *Nk;            
		float *sum_a;
		float *sum_b;
		float *d_fc_ball;
		//xyj��cpp������-ylj��F90������-���׶�Ӧ
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
		int sits_mode = 0;//ѡ��SITSģʽ
		int atom_numbers = 0;//��ϵ������ԭ����������һ��Ҫ��MD_INFORMATIONһ����������Ҫ
		int protein_atom_numbers = 0;//����SITS��ǿ��ԭ������ԭ������һ��Ҫ��[0-protein_atom_numbers)��ģ������ҪԤ������
		int max_neighbor_numbers = 800;//�����һ����neighbor_list�����һ��
		float pwwp_enhance_factor = 0.5; //�໥�����ܲ����ж��ټ�ǿ
        float fc_ball = 1.;//SITS��ǿ������
	}info;

	
	float *d_total_pp_atom_energy = NULL;//AA������
	float *d_total_ww_atom_energy = NULL;//BB������
	float h_total_pp_atom_energy = 0;//AA������
	float h_total_ww_atom_energy = 0;//BB������    
    
	float *d_protein_water_atom_energy = NULL;//���ڼ�¼ÿ��ԭ�ӵ�AB������
	float *d_total_protein_water_atom_energy = NULL;//AB������
    float h_total_protein_water_atom_energy = 0;//AB������
	VECTOR *protein_water_frc = NULL;//���ڼ�¼AB����ԭ�ӽ�����������
	ATOM_GROUP *d_nl_ppww = NULL;//A��A�ھ���B��B�ھ�
	ATOM_GROUP *d_nl_pwwp = NULL;//A��B�ھ���B��A�ھ�

	//��ʼ��
	void Initial_SITS(CONTROLLER *controller, int atom_numbers);
	//�����һ�����������������ȷ��һ���Ƿ���Ҫ��������
	void Prepare_For_Calculate_Force(int *need_atom_energy, int isPrintStep);
	
	//����ͳ������ѧԭ����ǿ����ֻ��Ҫ������ƥ��ͺã���˲���Ҫ�����е��໥���ý�����ǿ
	//20201217��Ŀǰ��Ҫ��ǿ������bond,angle,dihedral,LJ��PME_Direct,nb_14���������ֲ���ǿ����Ӱ������ȷ�ԣ�
	//�ڼ������SITS��ص�����֮ǰ�����һ��
	void Clear_SITS_Energy();
	


	//���������������м�Ϳ�������������Ҫ��ǿ������������ģ���ԭ���������㺯��
	//�ڼ������SITS������������������
	void Calculate_Total_SITS_Energy(float *d_atom_energy);
    
	//�ı�frc��ʹ����Ҫ��ǿ��frc��ѡ������ǿ�����ڹ��õ�frc������ⲽfrc��ǿ��Ҫ�ŵ��պü���������Ҫ��ǿ��frc�ĺ����·�����������ǿ��Ӧ����ǿ����ģ��frc
	//��ˣ���ϵ�����п��ܵ�frc��Ҫ�������ǿ��frc������˺��������㲻Ӧ��ǿ��frc
	void SITS_Enhanced_Force(int steps, VECTOR *frc);

	//����Ϣ��ӡ
	void Print();

	//����������ͷ�SITS�ڵĸ���ָ��
	void Clear_SITS();
};
#endif //SITS_CUH(sits.cuh)
