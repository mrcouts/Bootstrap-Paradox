#include "RK.h"
#include <iostream>

RK::RK(string method, vec (*f_)(double, vec)){
	caso = 1;
	this->f_ = f_;
	this->SetMethod(method);
}

RK::RK(string method, Acceleration *AC){
	caso = 2;
	this->AC = AC;
	this->SetMethod(method);
}

RK::RK(string method, GNR *gnr){
	caso = 3;
	this->gnr = gnr;
	this->SetMethod(method);
}

RK::RK(string method, Serial *R, ControlLaw *CL, int nh){
	caso = 4;
	this->R = R;
	this->CL = CL;
	this->nh = nh;
	this->SetMethod(method);
}

void RK::SetMethod(string method){
	if(method == "Heun"){
		N = 2;
		a__ = new vec[N-1];
		a__[0] = {1.0};
		b_ = {0.5, 0.5};
	}
	else if(method == "RK3"){
		N = 3;
		a__ = new vec[N-1];
		a__[0] = {0.5};
		a__[1] = {-1.0, 2.0};
		b_ = {1.0/6, 2.0/3, 1.0/6};
	}
	else if(method == "RK4"){
		N = 4;
		a__ = new vec[N-1];
		a__[0] = {0.5};
		a__[1] = {0, 0.5};
		a__[2] = {0, 0, 1.0};
		b_ = {1.0/6, 1.0/3, 1.0/3, 1.0/6};
	}
	else if(method == "RK6"){
		N = 7;
		a__ = new vec[N-1];
		a__[0] = {0.5};
		a__[1] = {2.0/9, 4.0/9};
		a__[2] = {7.0/36, 2.0/9, -1.0/12};
		a__[3] = {-35.0/144, -55.0/36, 35.0/48, 15.0/8};
		a__[4] = {-1.0/360, -11.0/36, -1.0/8, 0.5, 1.0/10};
		a__[5] = {-41.0/260, 22.0/13, 43.0/156, -118.0/39, 32.0/195, 80.0/39};
		b_ = {13.0/200, 0, 11.0/40, 11.0/40, 4.0/25, 4.0/25, 13.0/200};
	}
	else if(method == "RK8"){
		N = 13;
		a__ = new vec[N-1];
		a__[0] = {0.5555555555555555555555555555555555555555555555555555555555555555555555555555555555556e-1};
		a__[1] = {0.2083333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1, 
			      0.625e-1};
		a__[2] = {0.3125e-1, 0, 0.9375e-1};
		a__[3] = {0.3125, 0, -1.171875, 1.171875};
		a__[4] = {0.375e-1, 0, 0, 0.1875, 0.15};
		a__[5] = {0.4791013711111111111111111111111111111111111111111111111111111111111111111111111111111e-1, 
                  0,
                  0,
                  0.1122487127777777777777777777777777777777777777777777777777777777777777777777777777778,
                  0.2550567377777777777777777777777777777777777777777777777777777777777777777777777777778e-1,
                  0.1284682388888888888888888888888888888888888888888888888888888888888888888888888888889e-1};
        a__[6] = {0.1691798978729228118143110713603823606551492879543441068765183916901985069878116840258e-1, 
                  0,
                  0,
                  0.3878482784860431695265457441593733533707275558526027137301504136188413357699752956243,
                  0.3597736985150032789670088963477236800815873945874968299566918115320174598617642813621e-1, 
                  0.1969702142156660601567152560721498881281698021684239074332818272245101975612112344987, 
                 -0.1727138523405018387613929970023338455725710262752107148467532611655731300161441266618};
        a__[7] = {0.6909575335919230064856454898454767856104137244685570409618590205329379141801779261271e-1,  
                  0,
                  0,
                 -0.6342479767288541518828078749717385465426336031120747443029278238678259156558727343870, 
                 -0.1611975752246040803668769239818171234422237864808493434053355718725564004275765826294, 
                  0.1386503094588252554198669501330158019276654889495806914244800957316809692178564181463,  
                  0.9409286140357562697242396841302583434711811483145028697758469500622977999086075976730, 
                  0.2116363264819439818553721171319021047635363886083981439225363020792869599016568720713};
        a__[8] = {0.1835569968390453854898060235368803084497642516277329033567822939550366893470341427061, 
                  0,
                  0,
                 -2.468768084315592452744315759974107457777649753547732495587510208118948846864075671103,
                 -0.2912868878163004563880025728039519800543380294081727678722180378521228988619909525814,  
                 -0.2647302023311737568843979946594614325963055282676866851254669985184857900430792761417e-1, 
                  2.847838764192800449164518254216773770231582185011953158945518598604677368771006006947,
                  0.2813873314698497925394036418267117820980705455360140173438806414700945017562775967023,
                  0.1237448998633146576270302126636397203122013536069738523260934117931117648560568049432};
        a__[9] ={-1.215424817395888059160510525029662994880024988044112261492933435562347094075572196306, 
                  0,
                  0,
                  16.67260866594577243228041328856410774858907460078758981907659463146817850515953150318,
                  0.9157418284168179605957186504507426331593736334298114556675054711691705693180672015549, 
                 -6.056605804357470947554505543091634004081967083324413691952297768849489422976536037188, 
                 -16.00357359415617811184170641007882303068079304063907122528682690677051264091851070681,
                  14.84930308629766255754539189802663208272299893302745636963476380528935538206408294404,
                 -13.37157573528984931829304139618159579089195289469740214314819172008509426917249413996,
                  5.134182648179637933173253611658602898712494286162881495270691720760048063805245199662};
        a__[10] ={0.2588609164382642838157309322317577667296307766301063163257807024364127266506000943372,  
                  0,
                  0,
                 -4.774485785489205112310117509706042746829391853746564335432052648850273207666153414325,
                 -0.4350930137770325094407004118103177819323551661617974116300885410959889587870529265443, 
                 -3.049483332072241509560512866312031613982854911220735245095857885009959537455491093975, 
                  5.577920039936099117423676634464941858623588944531498679305832747426169795721583100328,
                  6.155831589861040097338689126688954481197754937462937466615262374558053507207577588692,
                 -5.062104586736938370077406433910391644990220712141673880668482097753119682588189892664,
                  2.193926173180679061274914290465806019788262707389033759251111926114017568440058142981,
                  0.1346279986593349415357262378873236613955852772571946513284934221746877884770684011715};
        a__[11] ={0.8224275996265074779631682047726665909572303617765850630130165537026064642659285212794, 
                  0,
                  0,
                 -11.65867325727766428397655303545841477547369082638864247596731917545919361630644756876, 
                 -0.7576221166909361958811161540882449653663757591954118634754438929149012424414834787423,  
                  0.7139735881595815279782692827650546753142248878566301594716125352788068216039494192152,  
                  12.07577498689005673956617044860067967095705800972197897084832370633043082304810801069,
                 -2.127659113920402656390820858969398635427927973275119750336406473203376323696576615278,
                  1.990166207048955418328071698344314152176173017359794982793519250552799420955298804232,
                 -0.2342864715440402926602946918568015314512425817004394700255034276765120670064339827205,
                  0.1758985777079422650731051058901448183145508638446243836782009233893397195776568900819, 
                  0};
		b_ = {0.4174749114153024622208592846850711513419360280744778485846131359556941056165592859127e-1,
			  0,
			  0,
			  0,
			  0,
			 -0.5545232861123930896152189465471671889359328156519226303661201143086248409346173220672e-1,
			  0.2393128072011800970467473542487569696603053593110659071772286813575015789696127063297,
			  0.7035106694034430230580464108897021513663799729422329988912282232809616439130702896567,
			 -0.7597596138144609298844876770850584076554230192157088121727371293320966312254783676226,
			  0.6605630309222863414613785948378206399404197125341442757648140724130090459250911959182,
			  0.1581874825101233355296148386006854439728567324034642678211136427219891617064005758134,
			 -0.2381095387528628044718635553056971935251390792174541593034967926060717257568905964800,
			  0.25};
	}
	else cout << "Método indisponível!" << endl;

	c_.zeros(N-1);
	for(int i = 0; i<N-1; i++)
		c_(i) = sum(a__[i]);
}

RK::~RK(){
	delete[] a__;
	b_.clear() ;
	c_.clear();

	t_.clear();
	y__.clear();
	dy__.clear();
	u__.clear();
	k__.clear();
}

void RK::Doit(double h, double tf, vec y0_){
	int nt = int(tf/h);
	t_.zeros(nt+1);
	for(int i = 0; i<nt+1; i++)
		t_(i) = i*h;
	y__.zeros(y0_.n_rows, 1, nt+1);
	dy__.zeros(y0_.n_rows, 1, nt+1);
	u__.zeros(y0_.n_rows/2, 1, nt+1);
	y__.slice(0) = y0_;
	k__.zeros(y0_.n_rows, 1, N);
	counter = 0;

	vec aux_; aux_ = zeros(y0_.n_rows);
	for(int i = 0; i<nt; i++){
		field<vec> F2(2);
		switch (caso){
		    case 1: k__.slice(0) = f_(t_(i),y__.slice(i)); break;
		    case 2:
		        F2 = AC->f2_(t_(i),y__.slice(i));
		        for(uint j = 0; j< F2(1).n_rows; j++) u__(j,0,i) = F2(1)(j);
		        //u__.slice(i) = F2(1); 
		        k__.slice(0) = F2(0);
		        break;
		    case 3: k__.slice(0) = gnr->g_(y__.slice(i)); break;
		    case 4:
		    	if(counter == 0){
		    		k__.slice(0) = R->f_(y__.slice(i), CL->Doit(t_(i), y__(span(0,R->dof-1),span(0,0),span(i,i)), dy__(span(0,R->dof-1),span(0,0),span(i,i)) ) );
		    		counter = nh - 1;
		    	}
		    	else{
		    		k__.slice(0) = R->f_(y__.slice(i), R->u_ );
		    		counter--;
		    	}
		    	for(uint j = 0; j< R->dof; j++) u__(j,0,i) = R->u_(j);
		    	break;
		}
		for(int j = 1; j<N; j++){
			aux_.zeros();
			for(int k = 0; k<j; k++)
				aux_ += a__[j-1](k)*k__.slice(k);
			switch (caso){
		        case 1: k__.slice(j) = f_(t_(i) + c_(j-1)*h, y__.slice(i) + h*aux_); break;
		        case 2: k__.slice(j) = AC->f_(t_(i) + c_(j-1)*h, y__.slice(i) + h*aux_); break;
		        case 3: k__.slice(j) = gnr->g_(y__.slice(i) + h*aux_); break;
		        case 4: k__.slice(j) = R->f_(y__.slice(i) + h*aux_, R->u_); break;
		    }
		aux_.zeros();
		for(int i = 0; i<N; i++)
			aux_ += b_(i)*k__.slice(i);
		y__.slice(i+1) = y__.slice(i) + h*aux_; 
		dy__.slice(i+1) = (y__.slice(i+1) - y__.slice(i))/h ; 
	    }
	}
	if(caso == 2){
		vec aux2_ = AC->f2_(t_(nt),y__.slice(nt))(1);
		for(uint j = 0; j< aux2_.n_rows; j++) u__(j,0,nt) = aux2_(j);
		//u__.slice(nt) = AC->f2_(t_(nt),y__.slice(nt))(1);
	}
	if(caso == 4){
		vec aux2_ = R->f_(y__.slice(nt), CL->Doit(t_(nt), y__(span(0,R->dof-1),span(0,0),span(nt,nt)), y__(span(R->dof,2*R->dof-1),span(0,0),span(nt,nt)) ) );
		for(uint j = 0; j< R->dof; j++) u__(j,0,nt) = R->u_(j);
	}
}