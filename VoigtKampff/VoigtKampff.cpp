#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include "VoigtKampff.h"
#include <cstdlib>
#include <cstdio>
const double PI = 3.14159265359;
const double SQRTLN2 = 0.832554611;
const double ISQRTPI = 0.564189584; 

const double DISTANCE_MAGIC_NUMBER = 1.9;
const double REFERNECE_NU = 1.0;
double CubicInterpolate(
   double y0,double y1,
   double y2,double y3,
   double mu)
{
   double a0,a1,a2,a3,mu2;

   mu2 = mu*mu;
   a0 = y3 - y2 - y0 + y1;
   a1 = y0 - y1 - a0;
   a2 = y2 - y0;
   a3 = y1;

   return(a0*mu*mu2+a1*mu2+a2*mu+a3);
};


VoigtKampff::VoigtKampff(double pGammaD,double pRes,double pLorenzCutoff,bool pNormalize){
	m_Npoints = 2*pLorenzCutoff/pRes + 1;
	middle_point = m_Npoints/2;
	m_res=pRes;
	m_gammaD = pGammaD;
	m_lorentz_cutoff = pLorenzCutoff;
	m_normalize = pNormalize;
	num_gammaL=0;


}


VoigtKampff::VoigtKampff(std::vector<double> pGammaL,double pGammaD,double pRes,double pLorenzCutoff,bool pNormalize) : VoigtKampff(pGammaD,pRes,pLorenzCutoff,pNormalize){
	int dummy;
	//num_gammaL=0;
	//m_Npoints = 2*pLorenzCutoff/pRes + 1;
	//middle_point = m_Npoints/2;
	//m_res=pRes;
	//m_gammaD = pGammaD;
	//m_lorentz_cutoff = pLorenzCutoff;
	
	std::cout<<"Initializing Voigt....."<<std::endl;
	std::cout<<"Npoints: "<<m_Npoints<<" middle: "<<middle_point<<"Humlecik points: "<<DISTANCE_MAGIC_NUMBER/m_res<<" range : "<< 0.0 << " - "<< double(m_Npoints)*m_res<<std::endl;
	m_gammaL = pGammaL;
	double x,y;
	double nu = 0.0;
	for(int i = 0; i < pGammaL.size(); i++){
		double gammaL = pGammaL[i];
		AddLorentzian(gammaL,&dummy);
	}

	std::cout<<"Done!!!"<<std::endl;



}

int VoigtKampff::CheckLorentzian(double gammaL){
			int gammaIdx=-1;
			std::map<float,int>::iterator i = m_gamma_mapping.find((float)gammaL);
			if(i == m_gamma_mapping.end()){
				gammaIdx = -1;
				//AddLorentzian(gammaL);
			}else{
				gammaIdx = i->second;
			}	
			return gammaIdx;
}


double VoigtKampff::ComputeVoigt(double dfreq,int gammaL,double nu){
	
	double dfreq_=fabs(dfreq);
	if(dfreq_ > m_lorentz_cutoff)
		return 0.0;
	int point = dfreq/m_res + middle_point;
	//printf("point:%d\n",point);
	//point = std::max(point,0c);
	//point = std::min(point,m_Npoints);
	
		
		//TODO: Use intepolation 
		/*int i1 = std::max(point-2,0);
		int i2 = std::max(point-1,0);
		int i3 =  std::min(point+2,m_Npoints);
		int i4 = std::min(point+2,m_Npoints);
		double p1=m_voigt_grid[gammaL][i1];
		double p2=m_voigt_grid[gammaL][i2];
		double p3=m_voigt_grid[gammaL][i3];
		double p4=m_voigt_grid[gammaL][i4];
		return CubicInterpolate(p1,p2,p3,p4,dfreq);///ComputeDoppler(dfreq,m_gammaD,nu);
		*/
	//double dist = 2.0*m_gammaD*nu/m_gammaL[gammaL];
	//if(dfreq_ < dist) return HumlicekTest(dfreq,m_gammaL[gammaL],nu);
	return m_voigt_grid.at(gammaL).at(point);
		//;
	//return 0.0;
	
}	

void VoigtKampff::AddLorentzian(double gammaL,int* index){
		double nu = -m_lorentz_cutoff;
		*index = m_voigt_grid.size();
		m_voigt_grid.push_back(std::vector<double>());	
		
		m_gamma_mapping[(float)gammaL]= num_gammaL;
		m_gammaL.push_back(gammaL);
		m_mag.push_back(0.0);
		for(int j = 0; j < m_Npoints; j++){	
			
			double x = SQRTLN2*fabs(nu)/(m_gammaD*REFERNECE_NU);
			double y = SQRTLN2*gammaL/(m_gammaD*REFERNECE_NU);
		
			m_voigt_grid.back().push_back(humlic(x,y)*SQRTLN2*ISQRTPI/(m_gammaD*REFERNECE_NU));
			m_mag.back() += m_voigt_grid.back().back();
			nu+=m_res;
		}
		std::cout<<"mag:"<< m_mag.back()*m_res<<std::endl;
		m_mag.back() *= m_res;
		num_gammaL++;
}


void VoigtKampff::DoVectorized(double* __restrict intens,const std::vector<double> & gammaVoigt,const double abscoef,const int start,const int end){
	//int count = end-start;
	
	for(int i = start; i < end; i++){
		intens[i]+=gammaVoigt[i]*abscoef;
	}


}

void VoigtKampff::ComputeVoigtVectorized(const double* __restrict freq,double* __restrict intens,const double abscoef,const int ib,const int ie,const double gammaL,const double start_nu,const double nu){
			int gammaIdx = CheckLorentzian(gammaL);
			//printf("%d\n",gammaIdx);


			ComputeVoigtVectorized(freq, intens,abscoef,ib,ie,gammaIdx,start_nu,nu);
}

void VoigtKampff::ComputeVoigtVectorized(const double* __restrict freq,double* __restrict intens,const double abscoef,const int ib,const int ie,const int gammaL,const double start_nu,const double nu){


	//int ib=ib_-1;
	int center_point = (nu-start_nu)/m_res;
	int middle_shift = center_point - middle_point;
	double gammaL_ = m_gammaL[gammaL];
	
	//int dist = DISTANCE_MAGIC_NUMBER/m_res;//(m_gammaD*nu/m_gammaL[gammaL])/m_res;
	int dist;
	if (gammaL_ == 0.0)
		dist = m_Npoints;
	else
		//dist = std::max(m_gammaD / m_res,2.0);
		//dist = std::max((m_gammaD*nu/m_gammaL[gammaL])/m_res,8.0);
		dist = std::max(2.0/m_res,3.0);
	int ib_rel = ib - middle_shift;
	int ie_rel = ie - middle_shift;
	
	double gammaD_ = (m_gammaD*nu);

	
	//dist = std::max(101,dist);
	//int start_dist 
	//printf("Middle calc %d %d\n",start_dist,end_dist);
	
	
	
	double dfreq;
	//DO PURE DOPPLER
	if (gammaL_ == 0.0) {
		//USe vectorized doppler
		DoDopplerVectorized(freq, intens, abscoef, ib, ie, gammaD_, nu);



	}//OTHERWISE WE VOIGT BABEEEEEEEEEEEEEEE
	else {
		int Npoints = ie - ib + 1;
		double mag = m_mag[gammaL];

		if(!m_normalize) mag = 1.0;
		
		int begin_humlicek;
		int end_humlicek;
		int begin_reflect;
		int end_reflect;
		int fudge_point;
		
		//Do the middle
		int start_dist = std::max(center_point - dist, ib);
		int end_dist = std::min(center_point + dist, ie);
		
		int left_start = ib_rel;
		int left_end = std::min(middle_point - dist, ie_rel);
		int right_start = std::max(middle_point + dist, ib_rel);
		int right_end = ie_rel;
		
		if((left_end - left_start) >= (right_start - right_end)){
			begin_humlicek =start_dist;
			end_humlicek=std::min(center_point+2,ie);
			fudge_point = start_dist;
			begin_reflect=end_humlicek;
			end_reflect = end_dist;
			
		}else{
			begin_humlicek =std::max(center_point-1,start_dist);
			end_humlicek=end_dist;
			fudge_point = end_humlicek;
			begin_reflect=start_dist;
			end_reflect =begin_humlicek;
		}
		
		int total_hum_points;
		//printf("Middle calc %d %d %d %d %d %d %d %12.6f\n",begin_humlicek,center_point,end_humlicek,begin_reflect,end_reflect,ib,ie,freq[center_point]);
		
		double fudge_factor=1.0f;
		/*if (ib < center_point) {
			if(left_end >= left_start && left_end >0){
			 	fudge_factor =  HumlicekTest(freq[ib - ib_rel+left_end-1] - nu, gammaL_, nu)/m_voigt_grid[gammaL][left_end-1];
			 }

		
		}else if (ie >= center_point){
			if(right_end >= right_start && right_start >0){
				fudge_factor =  HumlicekTest(freq[ib - ib_rel+right_start-1] - nu, gammaL_, nu)/m_voigt_grid[gammaL][right_start-1];
			}
	
		}*/
		
		//mag /= fudge_factor;
		
		double* temp_humlicek;
		if (start_dist < end_dist) {
			total_hum_points = end_dist - start_dist + 1;
			temp_humlicek = new double[total_hum_points];

			for (int i = start_dist; i < end_dist; i++) {
				dfreq = freq[i] - nu;
				double hum = HumlicekTest(dfreq, gammaL_, nu);
				temp_humlicek[i - start_dist] = hum;
				mag -= m_voigt_grid[gammaL][i -middle_shift]*m_res;
				mag += hum*m_res;
			}
			
			//fudge_factor = temp_humlicek[fudge_point - start_dist]/m_voigt_grid[gammaL][fudge_point -middle_shift];
			fudge_factor = 1.0;
			//This is so that the second one is vectorized
			if(!m_normalize) mag = 1.0;
			
			for(int i = start_dist; i < end_dist; i++){
				intens[i] += temp_humlicek[i-start_dist]*abscoef/mag;
				
				
			}
			
			delete[] temp_humlicek;

		}
		//printf("ib: %d ie: %d start:%d end:%d dist: %d center:%d nu: %12.6f\n",ib,ie,ib_rel,ie_rel,dist,center_point,nu);

		//If we have a calculation on the left
		if (ib < center_point) {
			DoVectorized(intens + ib - ib_rel, m_voigt_grid[gammaL], (abscoef/mag)*fudge_factor, left_start, left_end);
			//printf("Left calc %d %d\n",left_start,left_end);

		}
		if (ie >= center_point) {

			DoVectorized(intens + ib - ib_rel, m_voigt_grid[gammaL], (abscoef/mag)*fudge_factor, right_start, right_end);
			//printf("Right calc %d %d\n",right_start,right_end);

		}


	}
	


	//int dist_hum = dist;
	
	//for(int i = dist+1; i < dist
	
	
	
	//exit(0);
	
	
	
	



}

void VoigtKampff::DoDopplerVectorized(const double * __restrict freq, double * __restrict intens, const double abscoef, const int ib, const int ie, const double gammaD, const double nu)
{
	double fact = SQRTLN2 / gammaD;
	double x0 = fact*m_res*0.5;
	
	for (int i = ib; i < ie; i++) {
	//	double dfreq = freq[i] - nu;
		//double xp = fact*(dfreq)+x0;
		//double xm = fact*(dfreq)-x0;
	//	double de = erf(xp) - erf(xm);
		intens[i] += abscoef*0.5/m_res*(erf(fact*(freq[i] - nu) + x0)
			- erf(fact*(freq[i] - nu) - x0));
			//intens[i] += exp(i);
			/*intens[i] +=
				//DE
				( erf(   fact*(freq[i] - nu) + x0   )
				- erf(  fact*(freq[i] - nu) - x0  ) )

				*abscoef*0.5 /
				(freq[i] - nu);
		}*/
	}
}


double VoigtKampff::ComputeDoppler(double dfreq,double gammaD,double nu){
	
    	 //double dop = 0.03912;//(3155.87/nu)*4.523*(gammaD*nu);
         /*  
         double xp = SQRTLN2/dop*(freq+ m_res/2.0 - nu);
         double xm = SQRTLN2/dop*(freq- m_res/2.0 - nu);
             
         double de = erf(xp)-erf(xm);
        
         return 0.5/m_res*de;
	*/
	
	//double alpha = -1.0/(dop);
	//double de = exp(alpha*dfreq_*dfreq_);
	//double de = (dfreq*dfreq + (dop)*(dop))/(dop*dop*nu) ;
	//return de;
	return (-3.295e-05*dfreq*dfreq*dfreq*dfreq + 2.916e-08*dfreq*dfreq*dfreq 
		+ 0.2924*dfreq*dfreq -6.335e-06*dfreq + 0.02171);


}


	//void Initialize()
double VoigtKampff::ComputeVoigt(double dfreq,double gammaL){
	
	double dfreq_=fabs(dfreq);
	if(dfreq_ > m_lorentz_cutoff)
		return 0.0;
	int point = fabs(dfreq_)/m_res;
	//point = std::max(point,0);
	
	
	if(m_gamma_mapping.count(gammaL) > 0){
		int gamma_idx = m_gamma_mapping[gammaL];
		
		//TODO: Use intepolation 
		/*int i1 = std::max(point-1,0);
		int i2 = point;
		int i3 =  std::min(point+2,m_Npoints);
		int i4 = std::min(point+2,m_Npoints);
		double p1=m_voigt_grid[gamma_idx][i1];
		double p2=m_voigt_grid[gamma_idx][i2];
		double p3=m_voigt_grid[gamma_idx][i3];
		double p4=m_voigt_grid[gamma_idx][i4];
		return CubicInterpolate(p1,p2,p3,p4,dfreq);
		*/
		
		return m_voigt_grid[gamma_idx][point];
		//;
	}
	return 0.0;
	
}



double VoigtKampff::HumlicekTest(double dfreq,double gammaL,double nu){
	double x,y;
	double dfreq_=fabs(dfreq);
	double gammaD = m_gammaD*nu;
	
	x = SQRTLN2*dfreq_/(gammaD);
	y = SQRTLN2*gammaL/gammaD;
	return humlic(x,y)*SQRTLN2*ISQRTPI/gammaD;
}



