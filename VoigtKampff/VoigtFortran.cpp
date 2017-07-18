#include "VoigtKampff.h"
#include <vector>


std::vector<VoigtKampff> voigtkampff_vector;



extern "C" {
 void initalize_voigt_kampff_(double * pGammaD,double * pRes,double * pLorenzCutoff,int * index,int* norm);
 void check_lorentzian_(double * gammaL, int* index,int* gamma_idx);
 void compute_voigt_(double* freq,double* intens,const double* abscoef,int* ib,int* ie,int * gammaL,double * start_nu,double * nu, int* index);
 void add_lorentzian_(double* gammaL,int* index,int* gamma_idx);
}


void initalize_voigt_kampff_(double * pGammaD,double * pRes,double * pLorenzCutoff,int * index,int* norm){
	*index = voigtkampff_vector.size();
	
	voigtkampff_vector.push_back(VoigtKampff(*pGammaD,*pRes,*pLorenzCutoff,(bool)*norm));


}

void check_lorentzian_(double * gammaL, int* index,int* gamma_idx){
	
	float temp_gamma = *gammaL;
	
	*gamma_idx = voigtkampff_vector.at(*index).CheckLorentzian(temp_gamma);

}


void add_lorentzian_(double* gammaL,int* index,int* gamma_idx){
	voigtkampff_vector.at(*index).AddLorentzian(*gammaL,gamma_idx);
}


void compute_voigt_(double* freq,double* intens,const double* abscoef,int* ib,int* ie,int * gammaL,double * start_nu,double * nu, int* index){

	voigtkampff_vector.at(*index).ComputeVoigtVectorized(freq,intens,*abscoef,*ib,*ie,*gammaL,*start_nu,*nu);


}
