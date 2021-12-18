//#include "SITS_MD_force.cuh"
////this angle force calculate method is copied from AMBER16, may opt
//__global__ void SITS_Angle_Force_Energy(const int angle_numbers,
//	const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR scaler,
//	const int *atom_a, const int *atom_b, const int *atom_c,
//	const float *angle_k, const float *angle_theta0,
//	VECTOR *frc,float *energy)
//{
//	int angle_i = blockDim.x*blockIdx.x + threadIdx.x;
//	if (angle_i < angle_numbers)
//	{
//		int atom_i = atom_a[angle_i];
//		int atom_j = atom_b[angle_i];
//		int atom_k = atom_c[angle_i];
//
//		float da = angle_theta0[angle_i];
//		float df = angle_k[angle_i];
//
//		VECTOR ij_dr, kj_dr;
//		VECTOR fi, fk;
//		float rij_1, rkj_1;
//		float rik_1;
//		float cst;
//
//		float ant;
//		float dfw;
//		float cik;
//		float sth;
//		float cii;
//		float ckk;
//		ij_dr.x = ((int)(uint_crd[atom_i].uint_x - uint_crd[atom_j].uint_x)) * scaler.x;
//		ij_dr.y = ((int)(uint_crd[atom_i].uint_y - uint_crd[atom_j].uint_y)) * scaler.y;
//		ij_dr.z = ((int)(uint_crd[atom_i].uint_z - uint_crd[atom_j].uint_z)) * scaler.z;
//
//		kj_dr.x = ((int)(uint_crd[atom_k].uint_x - uint_crd[atom_j].uint_x)) * scaler.x;
//		kj_dr.y = ((int)(uint_crd[atom_k].uint_y - uint_crd[atom_j].uint_y)) * scaler.y;
//		kj_dr.z = ((int)(uint_crd[atom_k].uint_z - uint_crd[atom_j].uint_z)) * scaler.z;
//
//		rij_1 = ij_dr.x*ij_dr.x + ij_dr.y*ij_dr.y + ij_dr.z*ij_dr.z;
//		rij_1 = 1. / rij_1;
//
//		rkj_1 = kj_dr.x*kj_dr.x + kj_dr.y*kj_dr.y + kj_dr.z*kj_dr.z;
//		rkj_1 = 1. / rkj_1;
//		rik_1 = sqrtf(rij_1*rkj_1);
//
//		cst = ij_dr.x*kj_dr.x + ij_dr.y*kj_dr.y + ij_dr.z*kj_dr.z;
//		cst = cst * rik_1;
//
//		cst = fmaxf(-0.999999, fminf(cst, 0.999999));
//		/*if (cst < -0.999999)
//		{
//		cst = -0.999999;
//		}
//		else if (cst>0.999999)
//		{
//		cst = 0.999999;
//		}*/
//
//		ant = acosf(cst);
//
//		da = ant - da;
//		df = df*da;
//
//		dfw = -2.*df / sinf(ant);
//
//		cik = dfw * rik_1;
//		sth = dfw*cst;
//		cii = sth * rij_1;
//		ckk = sth * rkj_1;
//
//		fi.x = -cik*kj_dr.x + cii*ij_dr.x;
//		fi.y = -cik*kj_dr.y + cii*ij_dr.y;
//		fi.z = -cik*kj_dr.z + cii*ij_dr.z;
//
//		fk.x = -cik*ij_dr.x + ckk*kj_dr.x;
//		fk.y = -cik*ij_dr.y + ckk*kj_dr.y;
//		fk.z = -cik*ij_dr.z + ckk*kj_dr.z;
//
//		atomicAdd(&frc[atom_i].x, fi.x);
//		atomicAdd(&frc[atom_i].y, fi.y);
//		atomicAdd(&frc[atom_i].z, fi.z);
//
//		atomicAdd(&frc[atom_k].x, fk.x);
//		atomicAdd(&frc[atom_k].y, fk.y);
//		atomicAdd(&frc[atom_k].z, fk.z);
//
//		fi.x = -fi.x - fk.x;
//		fi.y = -fi.y - fk.y;
//		fi.z = -fi.z - fk.z;
//		atomicAdd(&frc[atom_j].x, fi.x);
//		atomicAdd(&frc[atom_j].y, fi.y);
//		atomicAdd(&frc[atom_j].z, fi.z);
//		atomicAdd(&energy[atom_i], df*da);//energy
//		//printf("%f\n", fi.x);
//	}
//}
//
//__global__ void SITS_Bond_Force_Energy(const int bond_numbers, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR scaler,
//	const int *atom_a, const int *atom_b, const float *bond_k, const float *bond_r0, VECTOR *frc,float *energy)
//{
//	int bond_i = blockDim.x*blockIdx.x + threadIdx.x;
//	if (bond_i < bond_numbers)
//	{
//		int atom_i = atom_a[bond_i];
//		int atom_j = atom_b[bond_i];
//		float tempk = bond_k[bond_i];
//		float tempr = bond_r0[bond_i];
//		float r1;
//		float ene;
//		VECTOR dr, f;
//		float tempf;
//
//		dr.x = ((int)(uint_crd[atom_i].uint_x - uint_crd[atom_j].uint_x)) * scaler.x;
//		dr.y = ((int)(uint_crd[atom_i].uint_y - uint_crd[atom_j].uint_y)) * scaler.y;
//		dr.z = ((int)(uint_crd[atom_i].uint_z - uint_crd[atom_j].uint_z)) * scaler.z;
//
//		r1 = norm3df(dr.x, dr.y, dr.z);
//		tempf = r1 - tempr;
//		ene = tempk*tempf*tempf;
//
//		tempf = tempf * tempk * 2 / r1;
//		f.x = tempf * dr.x;
//		f.y = tempf * dr.y;
//		f.z = tempf * dr.z;
//
//		atomicAdd(&frc[atom_i].x, -f.x);
//		atomicAdd(&frc[atom_i].y, -f.y);
//		atomicAdd(&frc[atom_i].z, -f.z);
//
//		atomicAdd(&frc[atom_j].x, f.x);
//		atomicAdd(&frc[atom_j].y, f.y);
//		atomicAdd(&frc[atom_j].z, f.z);
//
//		atomicAdd(&energy[atom_i], ene);//energy
//	}
//}
//
//__global__ void SITS_Dihedral_Force_Energy
//(const int dihedral_numbers, const UNSIGNED_INT_VECTOR *uint_crd, const VECTOR scaler,
//const int *atom_a, const int *atom_b, const int *atom_c, const int *atom_d, const int *ipn, const float *pk, const float *gamc, 
//const float *gams, const float *pn, VECTOR *frc,float *energy)
//{
//	int dihedral_i = blockDim.x*blockIdx.x + threadIdx.x;
//	if (dihedral_i < dihedral_numbers)
//	{
//		int atom_i = atom_a[dihedral_i];
//		int atom_j = atom_b[dihedral_i];
//		int atom_k = atom_c[dihedral_i];
//		int atom_l = atom_d[dihedral_i];
//
//		int temp_ipn = ipn[dihedral_i];
//
//		float temp_pk = pk[dihedral_i];
//		float temp_pn = pn[dihedral_i];
//		float temp_gamc = gamc[dihedral_i];
//		float temp_gams = gams[dihedral_i];
//
//		float tempf, xij, yij, zij, xkj, ykj, zkj, xkl, ykl, zkl, dx, dy, dz, gx, gy, gz;
//		float s, ap, cphi, sphi, cosnp, sinnp, ct, ct0, z1, z2, z12, z11, z22, df, dums;
//		float dc1, dc2, dc3, dc4, dc5, dc6, dr1, dr2, dr3, dr4, dr5, dr6, drx, dry, drz;
//		float fxi, fyi, fzi, fxj, fyj, fzj, fxk, fyk, fzk, fxl, fyl, fzl;
//
//
//		float ene;
//
//		xij = ((int)(uint_crd[atom_i].uint_x - uint_crd[atom_j].uint_x)) * scaler.x;
//		yij = ((int)(uint_crd[atom_i].uint_y - uint_crd[atom_j].uint_y)) * scaler.y;
//		zij = ((int)(uint_crd[atom_i].uint_z - uint_crd[atom_j].uint_z)) * scaler.z;
//		xkj = ((int)(uint_crd[atom_k].uint_x - uint_crd[atom_j].uint_x)) * scaler.x;
//		ykj = ((int)(uint_crd[atom_k].uint_y - uint_crd[atom_j].uint_y)) * scaler.y;
//		zkj = ((int)(uint_crd[atom_k].uint_z - uint_crd[atom_j].uint_z)) * scaler.z;
//		xkl = ((int)(uint_crd[atom_k].uint_x - uint_crd[atom_l].uint_x)) * scaler.x;
//		ykl = ((int)(uint_crd[atom_k].uint_y - uint_crd[atom_l].uint_y)) * scaler.y;
//		zkl = ((int)(uint_crd[atom_k].uint_z - uint_crd[atom_l].uint_z)) * scaler.z;
//
//		dx = yij * zkj - zij * ykj;
//		dy = zij * xkj - xij * zkj;
//		dz = xij * ykj - yij * xkj;
//		gx = zkj * ykl - ykj * zkl;
//		gy = xkj * zkl - zkj * xkl;
//		gz = ykj * xkl - xkj * ykl;
//
//		fxi = sqrt(dx * dx + dy * dy + dz * dz);
//		fyi = sqrt(gx * gx + gy * gy + gz * gz);
//		ct = dx * gx + dy * gy + dz * gz;
//
//		if (fxi >= 0.001)
//			z1 = 1.0 / fxi;
//		else
//			z1 = 0.0;
//
//		if (fyi >= 0.001)
//			z2 = 1.0 / fyi;
//		else
//			z2 = 0.0;
//
//		z12 = z1 * z2;
//
//		if (z12 != 0)
//			fzi = 1.0;
//		else
//			fzi = 0.0;
//
//		s = xkj * (dz * gy - dy * gz) + ykj *(dx * gz - dz * gx) + zkj * (dy *gx - dx * gy);
//
//		tempf = ct * z12;
//		if (tempf > 1.0)
//			tempf = 1.0;
//		if (tempf < -1.0)
//			tempf = -1.0;
//		tempf = acosf(tempf);
//		tempf = fabsf(tempf);
//		if (s < 0)
//			tempf = -tempf;
//
//
//		ap = 3.141592654 - tempf;
//
//
//
//		cphi = cosf(ap);
//		sphi = sinf(ap);
//
//		ct0 = temp_pn * ap;
//		cosnp = cosf(ct0);
//		sinnp = sinf(ct0);
//
//		tempf = 1e-18;
//		if (sphi < 0)
//			tempf = -tempf;
//
//		dums = sphi + tempf;
//
//		if (fabsf(dums) < 1e-6)
//			df = fzi * temp_gamc * (temp_pn - (((temp_ipn - 1) & 1) ^ 1) * temp_ipn + (((temp_ipn - 1) & 1) ^ 1) * temp_ipn * cphi);
//		else
//			df = fzi * temp_pn * (temp_gamc * sinnp - temp_gams * cosnp) / dums;
//
//		z11 = z1 * z1;
//		z22 = z2 * z2;
//		dc1 = -gx * z12 - cphi * dx * z11;
//		dc2 = -gy * z12 - cphi * dy * z11;
//		dc3 = -gz * z12 - cphi * dz * z11;
//		dc4 = dx * z12 + cphi * gx * z22;
//		dc5 = dy * z12 + cphi * gy * z22;
//		dc6 = dz * z12 + cphi * gz * z22;
//
//		dr1 = df * (dc3 * ykj - dc2 * zkj);
//		dr2 = df * (dc1 * zkj - dc3 * xkj);
//		dr3 = df * (dc2 * xkj - dc1 * ykj);
//		dr4 = df * (dc6 * ykj - dc5 * zkj);
//		dr5 = df * (dc4 * zkj - dc6 * xkj);
//		dr6 = df * (dc5 * xkj - dc4 * ykj);
//		drx = df * (-dc2 * zij + dc3 * yij + dc5 * zkl - dc6 * ykl);
//		dry = df * (dc1 * zij - dc3 * xij - dc4 * zkl + dc6 * xkl);
//		drz = df * (-dc1 * yij + dc2 * xij + dc4 * ykl - dc5 * xkl);
//
//		ene = (temp_pk + cosnp * temp_gamc + sinnp * temp_gams) * fzi;
//
//		fxi = -dr1;
//		fyi = -dr2;
//		fzi = -dr3;
//		fxj = dr1 - drx;
//		fyj = dr2 - dry;
//		fzj = dr3 - drz;
//		fxk = drx + dr4;
//		fyk = dry + dr5;
//		fzk = drz + dr6;
//		fxl = -dr4;
//		fyl = -dr5;
//		fzl = -dr6;
//
//
//		atomicAdd(&energy[atom_i], ene);//energy
//
//		atomicAdd(&frc[atom_i].x, fxi);
//		atomicAdd(&frc[atom_i].y, fyi);
//		atomicAdd(&frc[atom_i].z, fzi);
//		atomicAdd(&frc[atom_j].x, fxj);
//		atomicAdd(&frc[atom_j].y, fyj);
//		atomicAdd(&frc[atom_j].z, fzj);
//		atomicAdd(&frc[atom_k].x, fxk);
//		atomicAdd(&frc[atom_k].y, fyk);
//		atomicAdd(&frc[atom_k].z, fzk);
//		atomicAdd(&frc[atom_l].x, fxl);
//		atomicAdd(&frc[atom_l].y, fyl);
//		atomicAdd(&frc[atom_l].z, fzl);
//	}
//}
//
//__global__ void SITS_Dihedral_14_LJ_Force_With_Direct_CF_Energy(const int dihedral_14_numbers, const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR boxlength,
//	const int *a_14, const int *b_14, const float *lj_scale_factor, const float *cf_scale_factor, const float *LJ_type_A, const float *LJ_type_B, VECTOR *frc,float *energy)
//{
//	int dihedral_14_i = blockDim.x*blockIdx.x + threadIdx.x;
//	if (dihedral_14_i < dihedral_14_numbers)
//	{
//		int int_x;
//		int int_y;
//		int int_z;
//		UINT_VECTOR_LJ_TYPE r1, r2;
//		VECTOR dr;
//
//		float dr_1;
//		float dr_2;
//		float dr_4;
//		float dr_6;
//		float dr_8;
//		float dr_12;
//		float dr_14;
//		float frc_abs = 0.;
//		VECTOR temp_frc;
//
//		float lj_ene = 0.;
//		float cf_ene = 0.;
//
//		int x, y;
//		int atom_pair_LJ_type;
//		float A, B;
//
//		int atom_i = a_14[dihedral_14_i];
//		int atom_j = b_14[dihedral_14_i];
//
//		r1 = uint_crd[atom_i];
//		r2 = uint_crd[atom_j];
//		int_x = r2.uint_x - r1.uint_x;
//		int_y = r2.uint_y - r1.uint_y;
//		int_z = r2.uint_z - r1.uint_z;
//		dr.x = boxlength.x*int_x;
//		dr.y = boxlength.y*int_y;
//		dr.z = boxlength.z*int_z;
//
//		dr_1 = rnorm3df(dr.x, dr.y, dr.z);
//		dr_2 = dr_1*dr_1;
//		dr_4 = dr_2*dr_2;
//		dr_6 = dr_4*dr_2;
//		dr_12 = dr_6*dr_6;
//		dr_8 = dr_4*dr_4;
//		dr_14 = dr_12*dr_2;
//
//
//		//CF
//		float charge_i = r1.charge;
//		float charge_j = r2.charge;
//		float frc_cf_abs;
//
//		cf_ene = charge_i*charge_j*dr_1;
//		cf_ene *= cf_scale_factor[dihedral_14_i];
//		frc_cf_abs = -cf_ene * dr_2 ;
//
//
//		//LJ
//		y = (r2.LJ_type - r1.LJ_type);
//		x = y >> 31;
//		y = (y^x) - x;
//		x = r2.LJ_type + r1.LJ_type;
//		r2.LJ_type = (x + y) >> 1;
//		x = (x - y) >> 1;
//		atom_pair_LJ_type = (r2.LJ_type*(r2.LJ_type + 1) >> 1) + x;
//
//		A = LJ_type_A[atom_pair_LJ_type];
//		B = LJ_type_B[atom_pair_LJ_type];
//		frc_abs = -A * dr_14
//			+ B * dr_8;
//		lj_ene = 0.08333333*A * dr_12
//			- 0.1666666*B * dr_6;//LJ的A,B系数已经乘以12和6因此要反乘
//
//		lj_ene *= lj_scale_factor[dihedral_14_i];
//		frc_abs *= lj_scale_factor[dihedral_14_i];
//
//		frc_abs += frc_cf_abs;
//		temp_frc.x = frc_abs*dr.x;
//		temp_frc.y = frc_abs*dr.y;
//		temp_frc.z = frc_abs*dr.z;
//
//
//
//		atomicAdd(&frc[atom_j].x, -temp_frc.x);
//		atomicAdd(&frc[atom_j].y, -temp_frc.y);
//		atomicAdd(&frc[atom_j].z, -temp_frc.z);
//		atomicAdd(&frc[atom_i].x, temp_frc.x);
//		atomicAdd(&frc[atom_i].y, temp_frc.y);
//		atomicAdd(&frc[atom_i].z, temp_frc.z);
//
//		atomicAdd(&energy[atom_i], lj_ene+cf_ene);//energy
//	}
//}
//
//__global__ void SITS_LJ_Force_With_Direct_CF_Energy(
//	const int atom_numbers, const NEIGHBOR_LIST *nl, const int *atomnumbers_in_nl,
//	const UINT_VECTOR_LJ_TYPE *uint_crd, const VECTOR boxlength,
//	const float *LJ_type_A, const float *LJ_type_B, const float cutoff,
//	VECTOR *frc, const float pme_beta, const float sqrt_pi,float *energy)
//{
//	int atom_i = blockDim.x*blockIdx.x + threadIdx.x;
//	if (atom_i < atom_numbers)
//	{
//		int N = atomnumbers_in_nl[atom_i];
//		//int B = ceilf((float)N / blockDim.y);
//		int atom_j;
//		int int_x;
//		int int_y;
//		int int_z;
//		UINT_VECTOR_LJ_TYPE r1 = uint_crd[atom_i], r2;
//		VECTOR dr;
//		float dr_1;
//		float dr_2;
//		float dr_4;
//		float dr_6;
//		float dr_8;
//		float dr_12;
//		float dr_14;
//		float frc_abs = 0.;
//		VECTOR frc_lin;
//		VECTOR frc_record = { 0., 0., 0. };
//
//		//CF
//		float charge_i = r1.charge; //r1.charge;
//		float charge_j;
//
//		float dr_abs;
//		float beta_dr;
//		float frc_cf_abs;
//		
//		//
//
//		int x, y;
//		int atom_pair_LJ_type;
//		float A, B;
//
//		float lj_ene=0.;
//		float cf_ene=0.;
//		float ene_lin = 0.;
//
//		for (int j = threadIdx.y; j < N; j = j + blockDim.y)
//		{
//
//			atom_j = nl[atom_i].atom_serial[j];
//			r2 = uint_crd[atom_j];
//			//CF
//			charge_j = r2.charge;
//
//			int_x = r2.uint_x - r1.uint_x;
//			int_y = r2.uint_y - r1.uint_y;
//			int_z = r2.uint_z - r1.uint_z;
//			dr.x = boxlength.x*int_x;
//			dr.y = boxlength.y*int_y;
//			dr.z = boxlength.z*int_z;
//			dr_abs = norm3df(dr.x, dr.y, dr.z);
//			if (dr_abs < cutoff)
//			{
//				dr_1 = 1. / dr_abs;
//				dr_2 = dr_1*dr_1;
//				dr_4 = dr_2*dr_2;
//				dr_6 = dr_4*dr_2;
//				dr_8 = dr_4*dr_4;
//				dr_12 = dr_6*dr_6;
//				dr_14 = dr_12*dr_2;
//
//				y = (r2.LJ_type - r1.LJ_type);
//				x = y >> 31;
//				y = (y^x) - x;
//				x = r2.LJ_type + r1.LJ_type;
//				r2.LJ_type = (x + y) >> 1;
//				x = (x - y) >> 1;
//				atom_pair_LJ_type = (r2.LJ_type*(r2.LJ_type + 1) >> 1) + x;
//
//				A = LJ_type_A[atom_pair_LJ_type];
//				B = LJ_type_B[atom_pair_LJ_type];
//				frc_abs = -A * dr_14+ B * dr_8;
//				lj_ene = 0.08333333*A * dr_12- 0.1666666*B * dr_6;//LJ的A,B系数已经乘以12和6因此要反乘
//				ene_lin = ene_lin + lj_ene;
//
//				//CF
//				beta_dr = pme_beta*dr_abs;
//				frc_cf_abs = beta_dr *sqrt_pi * expf(-beta_dr*beta_dr) + erfcf(beta_dr);
//				frc_cf_abs = frc_cf_abs * dr_2 *dr_1;
//				frc_cf_abs = charge_i * charge_j*frc_cf_abs;
//
//				frc_abs = frc_abs - frc_cf_abs;
//
//				cf_ene = charge_i * charge_j * erfcf(beta_dr) *dr_1;
//				ene_lin = ene_lin + cf_ene;
//
//
//				frc_lin.x = frc_abs*dr.x;
//				frc_lin.y = frc_abs*dr.y;
//				frc_lin.z = frc_abs*dr.z;
//
//				frc_record.x = frc_record.x + frc_lin.x;
//				frc_record.y = frc_record.y + frc_lin.y;
//				frc_record.z = frc_record.z + frc_lin.z;
//
//				atomicAdd(&frc[atom_j].x, -frc_lin.x);
//				atomicAdd(&frc[atom_j].y, -frc_lin.y);
//				atomicAdd(&frc[atom_j].z, -frc_lin.z);
//			}
//
//		}//atom_j cycle
//		atomicAdd(&frc[atom_i].x, frc_record.x);
//		atomicAdd(&frc[atom_i].y, frc_record.y);
//		atomicAdd(&frc[atom_i].z, frc_record.z);
//
//		atomicAdd(&energy[atom_i], ene_lin);
//	}
//}
//
//__global__ void SITS_For_Enhanced_Force_Calculate_NkExpBetakU_1
//(const int k_numbers, const float *beta_k, const float *nk,float *nkexpbetaku,
//const float ene)
//{
//	float lin = beta_k[k_numbers-1];
//	for (int i = threadIdx.x; i < k_numbers; i = i + blockDim.x)
//	{
//		nkexpbetaku[i] = nk[i] * expf(-(beta_k[i] - lin) * ene);
//		//printf("%f %f\n", beta_k[i], nkexpbetaku[i]);
//	}
//}
//__global__ void SITS_For_Enhanced_Force_Calculate_NkExpBetakU_2
//(const int k_numbers, const float *beta_k, const float *nk, float *nkexpbetaku,
//const float ene)
//{
//	float lin = beta_k[0];
//	for (int i = threadIdx.x; i < k_numbers; i = i + blockDim.x)
//	{
//		nkexpbetaku[i] = nk[i] * expf(-(beta_k[i] - lin) * ene);
//		//printf("%f %f\n", beta_k[i], nkexpbetaku[i]);
//	}
//}
//__global__ void SITS_For_Enhanced_Force_Sum_Of_Above
//(const int k_numbers, const float *nkexpbetaku, const float *beta_k, float *sum_of_above)
//{
//	if (threadIdx.x == 0)
//	{
//		sum_of_above[0] = 0.;
//	}
//	__syncthreads();
//	float lin = 0.;
//	for (int i = threadIdx.x; i < k_numbers; i = i + blockDim.x)
//	{
//		lin = lin + beta_k[i]*nkexpbetaku[i];
//	}
//	atomicAdd(sum_of_above, lin);
//}
//__global__ void SITS_For_Enhanced_Force_Sum_Of_NkExpBetakU
//(const int k_numbers,const float *nkexpbetaku, float *sum_of_below)
//{
//	if (threadIdx.x == 0)
//	{
//		sum_of_below[0] = 0.;
//	}
//	__syncthreads();
//	float lin = 0.;
//	for (int i = threadIdx.x; i < k_numbers; i = i + blockDim.x)
//	{
//		lin = lin + nkexpbetaku[i];
//		//printf("%f\n", nkexpbetaku[i]);
//	}
//	atomicAdd(sum_of_below, lin);
//}
//__global__ void SITS_For_Enhanced_Force_Protein
//(const int protein_numbers, VECTOR *md_frc,const VECTOR *pw_frc,const float factor)
//{
//	float lin = 0.5*(factor + 1.);
//	for (int i = threadIdx.x; i < protein_numbers; i = i + blockDim.x)
//	{
//		md_frc[i].x = factor*(md_frc[i].x) + lin*pw_frc[i].x;
//		md_frc[i].y = factor*(md_frc[i].y) + lin*pw_frc[i].y;
//		md_frc[i].z = factor*(md_frc[i].z) + lin*pw_frc[i].z;
//	}
//}
//__global__ void SITS_For_Enhanced_Force_Water
//(const int protein_numbers,const int atom_numbers, VECTOR *md_frc, const VECTOR *pw_frc, const float factor)
//{
//	float lin = 0.5*(factor + 1.);
//	//printf("%f\n", factor);
//	for (int i = threadIdx.x + protein_numbers; i < atom_numbers; i = i + blockDim.x)
//	{
//		md_frc[i].x = md_frc[i].x + lin*pw_frc[i].x;
//		md_frc[i].y = md_frc[i].y + lin*pw_frc[i].y;
//		md_frc[i].z = md_frc[i].z + lin*pw_frc[i].z;
//	}
//}
//__global__ void SITS_Enhanced_Force
//(const int atom_numbers, const int protein_atom_numbers,
//VECTOR *md_frc, const VECTOR *pw_frc,
//const float *pp_ene, const float *pw_ene,
//const int k_numbers,float *nkexpbetaku,
//const float *beta_k, const float *n_k,
//float *sum_a,float *sum_b,float *factor,
//const float beta_0, const float pe_a, const float pe_b, const float fb_bias)
//{
//	float ene = pp_ene[0] + 0.5*pw_ene[0];
//	ene = pe_a * ene + pe_b;
//	//if (protein_atom_numbers != 0)//可能可以不用判断
//	//{
//		ene = (float)ene;
//		//printf("ene %f\n", ene);
//	//}
//
//		if (ene > 0)
//		{
//			SITS_For_Enhanced_Force_Calculate_NkExpBetakU_1 << <1, 64 >> >
//				(k_numbers, beta_k, n_k, nkexpbetaku, ene);
//		}
//		else
//		{
//			SITS_For_Enhanced_Force_Calculate_NkExpBetakU_2 << <1, 64 >> >
//				(k_numbers, beta_k, n_k, nkexpbetaku, ene);
//		}
//
//	SITS_For_Enhanced_Force_Sum_Of_NkExpBetakU << <1, 128 >> >
//			(k_numbers, nkexpbetaku, sum_b);
//
//	SITS_For_Enhanced_Force_Sum_Of_Above << <1, 128 >> >
//		(k_numbers, nkexpbetaku, beta_k, sum_a);
//
//
//	//这段代码如果不做0判断，会有未知bug，sum_b[0]本不会为0，但第一帧有可能为0；
//	factor[0] = sum_a[0] / sum_b[0] / beta_0;
//	if (!isinf(factor[0]) && !isnan(factor[0]) && (factor[0] > 0.2 * factor[1]) && (factor[0] < 5 * factor[1]) )
//	{
//		factor[1] = factor[0];
//	}
//	else
//	{
//		factor[0] = factor[1];
//	}
//	float fc = factor[0] + fb_bias;// fminf(fmaxf(factor[0], 0.6), 1.1);//factor[0];
//	/*if (steps > 1000)
//		fc *= 1.2;*/
//	
//	
////	printf("factor %e sum0 %e %e ene %f lfactor %e\n", fc, sum_a[0], sum_b[0], ene, factor[1]);
//	__syncthreads();
//	
//
//	//line
//	//fc = (ene - 20.) / 80./2. + 0.2;
//
//	SITS_For_Enhanced_Force_Protein << <1, 128 >> >
//		(protein_atom_numbers, md_frc, pw_frc, fc);
//	SITS_For_Enhanced_Force_Water << <1, 128 >> >
//		(protein_atom_numbers, atom_numbers, md_frc, pw_frc, fc);
//
//}