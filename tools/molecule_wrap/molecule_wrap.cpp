/*
	由夏义杰创建，用以将轨迹中的分子的坐标不会被周期性映射到两边（wrap）。
	使用方法：编译（gcc或使用codeblocks等简单编译）后，执行程序假设名字为molecule_wrap。无错误检查
	./molecule_wrap parm文件名 原本的轨迹文件名  盒子信息文件名   wrap以后的轨迹文件名  [每隔多少步打印处理结果]

*/


#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string.h>

struct VECTOR
{
    float x;
    float y;
    float z;
};



typedef std::vector<int> LINKAGE;
LINKAGE *links;
VECTOR *crds;
int bond_numbers;
int atom_numbers;
VECTOR boxlength;
VECTOR half_boxlength;
VECTOR boxangle;

void get_linkage(const char *file_name)
{
	FILE *parm = fopen(file_name,"r");

	printf("START READING BOND INFORMATION FROM AMBERFILE:\n");

	char temps[1000];
	int i, tempi, bond_with_hydrogen;
    int tempa, tempb;

	while (1)
	{
		if (!fgets(temps, 1000, parm))
			break;

		if (strcmp(temps, "%FLAG POINTERS                                                                  \n") == 0
			|| strcmp(temps, "%FORMAT(10I8)\n") == 0)
		{
			fgets(temps, 1000, parm);

            fscanf(parm, "%d", &atom_numbers);
            printf("	atom_numbers is %d\n", atom_numbers);
            fscanf(parm, "%d", &tempi);
			fscanf(parm, "%d", &bond_with_hydrogen);
			fscanf(parm, "%d", &bond_numbers);
			bond_numbers += bond_with_hydrogen;
			printf("	bond_numbers is %d\n", bond_numbers);
            links = new LINKAGE[atom_numbers];
            crds = new VECTOR[atom_numbers];
		}


		if (strcmp(temps, "%FLAG BONDS_INC_HYDROGEN                                                        \n") == 0)
		{
			printf("	reading bond_with_hydrogen %d\n", bond_with_hydrogen);
			fgets(temps, 1000, parm);
			for (i = 0; i < bond_with_hydrogen; i++)
			{
				fscanf(parm, "%d\n", &tempa);
				fscanf(parm, "%d\n", &tempb);
				fscanf(parm, "%d\n", &tempi);
				tempa /= 3;
				tempb /= 3;
				tempi -= 1;
                links[tempa].push_back(tempb);
                links[tempb].push_back(tempa);

			}
		}

		if (strcmp(temps, "%FLAG BONDS_WITHOUT_HYDROGEN                                                    \n") == 0)
		{
			printf("	reading bond_without_hydrogen %d\n", bond_numbers - bond_with_hydrogen);
			fgets(temps, 1000, parm);
			for (i = bond_with_hydrogen; i < bond_numbers; i++)
			{
				fscanf(parm, "%d\n", &tempa);
				fscanf(parm, "%d\n", &tempb);
				fscanf(parm, "%d\n", &tempi);
				tempa /= 3;
				tempb /= 3;
				tempi -= 1;
                links[tempa].push_back(tempb);
                links[tempb].push_back(tempa);
			}
		}
	}
	printf("END(READING BOND INFORMATION FROM AMBERFILE)\n");
	fclose(parm);
}

void printf(LINKAGE *link, int number)
{
    for (int i = 0; i<number; i++)
    {
        printf("%d: ",i);
        for (int j = 0; j < link[i].size(); j++)
        {
            printf("%d ",link[i][j]);
        }
        printf("\n");
    }
}

std::vector<LINKAGE> molecules;

void get_molecule()
{
    bool *unsearched = new bool[atom_numbers];
    for (int i = 0; i< atom_numbers; i++)
        unsearched[i] = 1;
    LINKAGE queue;

    for (int i = 0; i< atom_numbers; i++)
    {
        if (unsearched[i])
        {
            LINKAGE mol;
            mol.push_back(i);
            unsearched[i] = 0;
            queue.insert(queue.end(),links[i].begin(),links[i].end());
            while (!queue.empty())
            {
                int temp = queue.back();
                queue.pop_back();
                if (unsearched[temp])
                {
                    queue.insert(queue.end(),links[temp].begin(),links[temp].end());
                    mol.push_back(temp);
                    unsearched[temp] = 0;
                }
            }
            molecules.push_back(mol);
        }
    }
}

void crd_wrap()
{
    LINKAGE mol;
    VECTOR crd_i;
    VECTOR crd_0;
	int atom_i;
    float dx, dy, dz;
    for (int i = 0; i<molecules.size(); i++)
    {
        mol = molecules[i];
        crd_0 = crds[mol[0]];
        //printf(&mol,1);
        for (int j = 1; j < mol.size(); j++)
        {
			atom_i = mol[j];
            crd_i = crds[atom_i];
            dx = crd_i.x - crd_0.x;
            dy = crd_i.y - crd_0.y;
            dz = crd_i.z - crd_0.z;
            if (dx > half_boxlength.x)
                crds[atom_i].x -= boxlength.x;
            else if (dx < -half_boxlength.x)
                crds[atom_i].x += boxlength.x;
            if (dy > half_boxlength.y)
                crds[atom_i].y -= boxlength.y;
            else if (dy < -half_boxlength.y)
                crds[atom_i].y += boxlength.y;
            if (dz > half_boxlength.z)
                crds[atom_i].z -= boxlength.z;
            else if (dz < -half_boxlength.z)
                crds[atom_i].z += boxlength.z;
        }
    }
}

int main(int argc, const char *argv[])
{

	int ntwx = 1000;
    char parm[256];
    char input[256];
    char mdbox[256];
    char output[256];
	if (argc <= 1)
	{
		printf("\n\tusage: \n\t\tmolecule_wrap -p parm_file -i input_dat_traj -box mdbox -o output_dat_traj\n");
		printf("\n\tdescription:\n\t\tFix the dat file to wrap the coordinates of atoms in the same molecule\n");
		return 0;
	}
    for (int i = 0; i<argc;i++)
    {
        if (strcmp(argv[i],"-h")==0)
        {
            printf("\n\tusage: \n\t\tmolecule_wrap -p parm_file -i input_dat_traj -box mdbox -o output_dat_traj\n");
			printf("\n\tdescription:\n\t\tFix the dat file to wrap the coordinates of atoms in the same molecule\n");
			return 0;
		}
        else if (strcmp(argv[i],"-p")==0)
        {
            i++;
            strcpy(parm,argv[i]);
        }
		else if (strcmp(argv[i],"-box")==0)
        {
            i++;
            strcpy(mdbox,argv[i]);
        }
		else if (strcmp(argv[i],"-i")==0)
        {
            i++;
            strcpy(input,argv[i]);
        }
		else if (strcmp(argv[i],"-o")==0)
        {
            i++;
            strcpy(output,argv[i]);
        }
    }
	

    get_linkage(parm);
    get_molecule();

    FILE *fr = fopen(input,"rb");
    FILE *fb = fopen(mdbox,"r");
    FILE *fw = fopen(output,"wb");

    int frame = 0;
    while (1)
    {
        frame += 1;
        size_t read_size =  fread(crds,sizeof(VECTOR),atom_numbers,fr);
        if (read_size==0)
            break;
        fscanf(fb,"%f %f %f %f %f %f",&boxlength.x,&boxlength.y,&boxlength.z,&boxangle.x,&boxangle.y,&boxangle.z);

        half_boxlength.x = 0.5 * boxlength.x;
        half_boxlength.y = 0.5 * boxlength.y;
        half_boxlength.z = 0.5 * boxlength.z;
        crd_wrap();
        fwrite(crds,sizeof(VECTOR),atom_numbers,fw);
        if (frame % ntwx == 0)
        {
            printf("Processed: %d\n",frame);
        }
    }
    fclose(fw);
    fclose(fr);


    return 0;
}
