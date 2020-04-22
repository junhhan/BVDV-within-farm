#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//---------------------------------------------------------------------------------------------------
// Structure Setup
typedef struct animal_node {
	int demo_group;		// 0: Calves_BH, 1: Calves_MA, 2: Young Heifers, 3: Breeding Heifers, 4: Mixed Aged cows, 5: Heifers/Stores, 6: Steers/bulls
	int age;			// Age (day old) of cattle.
	int sex;			// 0: Female, 1: Male.
	double weight;		// Weight of cattle
	int cull_code;		// Reason for culling. 0: Natural mortality (default), 1: Test & Cull, 2: Meat work, 3: Trade (MA only)
	int days_to_cull;	// Days left for culling due to empty/test positiveness.
	int is_pregnant;	// 0: No, 1: Yes.
	int is_aied;		// Conception via artificial insemination 0: No, 1: Yes.
	int calf_fate;		// Infection status of new-born calf
	int parity;			// Parity
	int days_to_calv;	// Gestation period of 281 days.
	int days_to_estrous;// Estrous return every U(18, 24) days.
	int is_milking;		// 0: No milking, 1: Milking.
	int inf_stat;		// 0: Susceptible, 1: TI, 2: Recovered, 3: iMmuned, 4: PI, 5: Vaccinated.
	int days_left;		// Days left; for the depletion of MAb U(112, 190), or to get recovered U(10, 20).
	int gest_at_inf;	// Gestation day at infection
	int is_vaccinated;	// 0: No, >= 1: Vaccination count.
	int vac_left;		// Days left for the conversion to susceptible.
	struct animal_node *dam_node;
	struct animal_node *next_node;
} animal_node;


//---------------------------------------------------------------------------------------------------
// VARIABLES
const int yr_burn_demo= 5;
const int yr_burn_bvdv= 0;
const int yr_sim= 5;
const int n_farms= 1;


//---------------------------------------------------------------------------------------------------
// FUNCTIONS
void _Update_date(int year[1], int day[1]);
int _Random_binomial(int n, double p);
double _Random_normal(double mu, double sd);
double _Random_uniform(double a, double b);
animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list);
void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);
void _Create_animal_list(int **farm_info, int **group_size, animal_node *animal_list[]);
int _Fate_of_calf(int inf_stat, int gest_at_inf);
void _Add_new_calf(int year[1], int **n_pi, int **group_size, int farm_id, animal_node *dam_node, animal_node *animal_list[]);
int _Vaccination(int inf_stat, int is_vaccinated);
int _ELISA(int inf_stat);
void _Test_cull(int farm_id, int **n_test, int pooled, int year[1], animal_node *current_node);
void _Update_demo(int year[1], int day[1], int **n_pi, int **farm_info, int **group_size, int **n_replace, int **n_purchase,
	int **n_vaccine, int **n_test, int *ag_tested, int *follow_up, int *SCR, int *VAC, int *DF, animal_node *animal_list[]);
void _Update_cull(int year[1], int **group_size, double **wgt_beef, int **n_surplus, animal_node *animal_list[]);
int _Abortion(int inf_stat, int days_to_calv);
void _Seed_infection(int year[1], int **group_size, int *n_troj, animal_node *animal_list[]);
void _Chage_seropositive(double *prop, int **group_size, animal_node *animal_list[]);
void _Update_infection(int year[1], int **farm_info, int **group_size, int **n_replace, double *betas, double *p_pi_neighbour, double *lambda, int **n_ti, int *DF, animal_node *animal_list[]);


//---------------------------------------------------------------------------------------------------
// MAIN FUNCTION
void outbreak_beef(int *seed, int *n_ma, double *betas, double *rho, double *p_pi_neighbour, double *lambda,
					int *SCR, int *VAC, int *DF,
					int *n_t, int *n_p, double *w_b, double *sup, double *pur, 
					double *c_test, double *c_vac, double *c_df) {

	// Take the variables
	srand((*seed)); // Seed for rand().
	//srand(time(NULL) + (*seed)); // Seed for rand().
	int *n_troj;
	n_troj= calloc(n_farms, sizeof(int));
	
	// Set the date variables
	int day[1], year[1];
	day[0]= 242; // One day before the PSC as we will update the date as the model iteration starts
	year[0]= -1;

	// Preparation	
	int i;
	int n_year= yr_burn_demo + yr_burn_bvdv + yr_sim;
	int d_burn_demo= 365 * yr_burn_demo;
	int d_burn_bvdv= 365 * yr_burn_bvdv;
	int d_sim= 365 * yr_sim;

	// Array to store farm management information	
	int **farm_info;
	farm_info= calloc(n_farms, sizeof(int *));
	// Array to record the number of replacement heifers to be selected for each farm
	int **n_replace;
	n_replace= calloc(n_farms, sizeof(int *));
	// Array to record group size for each farm
	int **group_size;
	group_size= calloc(n_farms, sizeof(int *));
	// Arrays to record production variables
	double **wgt_beef; // Produced live-weigth per production season
	wgt_beef= calloc(n_farms, sizeof(double *));
	int **n_surplus; 
	n_surplus= calloc(n_farms, sizeof(int *));
	int **n_purchase; 
	n_purchase= calloc(n_farms, sizeof(int *));
	// Arrays to record intervention variables
	int **n_vaccine; // Number of vaccinated cattle
	n_vaccine= calloc(n_farms, sizeof(int *));
	int **n_test; // Number of Ag ELISA tested cattle
	n_test= calloc(n_farms, sizeof(int *));
	// Indicator for follow-up test of the pooled test
	int *ag_tested;
	ag_tested= calloc(n_farms, sizeof(int));
	int *follow_up;
	follow_up= calloc(n_farms, sizeof(int));
	// Array to record disease variable
	int **n_pi; // Number of BVDV persistently infected animals
	n_pi= calloc(n_farms, sizeof(int *));
	int **n_ti; // Number of BVDV transiently infected animals
	n_ti= calloc(n_farms, sizeof(int *));


	for (i= 0; i< n_farms; i++) {
		farm_info[i]= calloc(5, sizeof(int));

		n_replace[i]= calloc(2, sizeof(int));
		group_size[i]= calloc(42, sizeof(int)); // For each farm, up to 42 (7 demo group X 6 infection status) demographic compartments can exist.
												//      calves_BH, calves_MA, young heifers, breeding heifers, mixed-age cows, heifers/stores, steers/bulls
												// S:   0,         1,         2,             3,                4,              5,              6      
												// T:   7,         8,         9,             10,               11,             12,             13      
												// R:   14,        15,        16,            17,               18,             19,             20
												// P:   21,        22,        23,            24,               25,             26,             27
												// M:   28,        29,        30,            31,               32,             33,             34
												// VI:  35,        36,        37,            38,               39,             40,             40
		wgt_beef[i]= calloc(n_year, sizeof(double));
		n_surplus[i]= calloc(n_year, sizeof(double));
		n_purchase[i]= calloc(n_year, sizeof(double));
		n_vaccine[i]= calloc(n_year, sizeof(int));
		n_test[i]= calloc(n_year, sizeof(int));
		n_pi[i]= calloc(n_year, sizeof(int));
		n_ti[i]= calloc(n_year, sizeof(int));

		n_troj[i]= (int)round((double)(*n_ma) * (*rho));

		farm_info[i][0]= (*n_ma); 	// Cow size (Median cow size based on the Dairy NZ)
		farm_info[i][1]= 328; 	// DOY of planned start of mating (PSM)
		farm_info[i][2]= 90; 	// DOY of weaning
		farm_info[i][3]= 59; 	// DOY of mixing
		farm_info[i][4]= 83; 	// DOY of pregnancy scan

		n_replace[i][0]= (int)round((double)farm_info[i][0] * 0.22); // Number of female calves that will be selected for replacement heifers at the initial demographic setting
		n_replace[i][1]= 0; // Target herd size of MA cows
	}
	//printf("Variables declared and initialised.\n");


	animal_node *animal_list[n_farms]; // Declare animal lists
	_Create_animal_list(farm_info, group_size, animal_list);
	//printf("Initial farm demographic set.\n");


	// Burn-in period to stabilise the demography
	for (i= 0; i< d_burn_demo; i++) {
		_Update_date(year, day);
		_Update_demo(year, day, n_pi, farm_info, group_size, n_replace, n_purchase, n_vaccine, n_test, ag_tested, follow_up, SCR, VAC, DF, animal_list);
		_Update_cull(year, group_size, wgt_beef, n_surplus, animal_list);
	} // Year 5

	// Burn-in period to stabilise BVDV
	for (i= 0; i< d_burn_bvdv; i++) {
		_Update_date(year, day);
		_Update_demo(year, day, n_pi, farm_info, group_size, n_replace, n_purchase, n_vaccine, n_test, ag_tested, follow_up, SCR, VAC, DF, animal_list);
		_Update_cull(year, group_size, wgt_beef, n_surplus, animal_list);
		_Update_infection(year, farm_info, group_size, n_replace, betas, p_pi_neighbour, lambda, n_ti, DF, animal_list);
	} // Year 15
	//_Chage_seropositive(prop, group_size, animal_list);

	// Simulation starts
	for (i= 0; i< d_sim; i++) {
		_Update_date(year, day);
		_Update_demo(year, day, n_pi, farm_info, group_size, n_replace, n_purchase, n_vaccine, n_test, ag_tested, follow_up, SCR, VAC, DF, animal_list);
		// BVDV incursion
		if (year[0] == yr_burn_demo + yr_burn_bvdv) {
			_Seed_infection(year, group_size, n_troj, animal_list);
		}
		_Update_infection(year, farm_info, group_size, n_replace, betas, p_pi_neighbour, lambda, n_ti, DF, animal_list);
		_Update_cull(year, group_size, wgt_beef, n_surplus, animal_list);
	} // Year 20


	// Update outcome variables
	double df_cost;
	df_cost= pow((3.2 * (double)(*n_ma) * 10000 / 2.84)  / acos(-1.0), 0.5) * 2.0 * acos(-1.0);
	
	for (i= 0; i< n_year; i++) {
		n_t[i]= n_ti[0][i];
		n_p[i]= n_pi[0][i];
		w_b[i]= wgt_beef[0][i] * 5.76;
		sup[i]= (double)n_surplus[0][i] * 1700.0;
		pur[i]= (double)n_purchase[0][i] * 43.15; // Balance between the price of buying replacement heifer (1300) and selling empty cow (1256.85)
		c_test[i]= (double)n_test[0][i] * 17.09;
		c_vac[i]= (double)n_vaccine[0][i] * 6.74;
		if ((*DF) == 1) {
			if (i == yr_burn_demo + yr_burn_bvdv) {
				c_df[i]= df_cost * 1.8;
			} else if (i > yr_burn_demo + yr_burn_bvdv) {
				c_df[i]= df_cost * 1.8 * 0.1;
			} else {
				c_df[i]= 0.0;
			}
		}
	}
	
	// Free the memory
	for (i= 0; i< n_farms; i++) {
		animal_node *eraser; // Declare a pointer to free the animal list
		eraser= animal_list[i];
		while(animal_list[i]) {
			eraser= animal_list[i];
			animal_list[i]= eraser->next_node;
			free(eraser);
		}
		free(farm_info[i]);
		free(n_replace[i]);
		free(group_size[i]);
		free(wgt_beef[i]);
		free(n_surplus[i]);
		free(n_purchase[i]);
		free(n_test[i]);
		free(n_vaccine[i]);
		free(n_pi[i]);
		free(n_ti[i]);
	}
	free(farm_info);
	free(n_replace);
	free(group_size);
	free(wgt_beef);
	free(n_surplus);
	free(n_purchase);
	free(n_vaccine);
	free(n_test);
	free(n_pi);
	free(n_ti);
	free(ag_tested);
	free(follow_up);
	free(n_troj);
}

//---------------------------------------------------------------------------------------------------
// FUNCTIONS IN DETAILS
// Update Date: It updates the calendar year and week.
void _Update_date(int year[1], int day[1]) {
	day[0] +=1;
	if (day[0] == 365) {day[0] %=365;}
	if (day[0] == 243) {year[0] +=1;} // Day 262 is PSC.
}

// Number of success from binomial distribution: It retunrs a number of success according to a probability "p".
int _Random_binomial(int n, double p) {
    int index;
	int n_pos= 0; 
	int i= 0;
	if (p > 1) {p= 1.0;}
	if (p < 0) {p= 0.0;}
	
	double x;

	if (p == 0.0) {
		return 0;
	} else if (p == 1.0) {
		return n;
	} else {
	    do {
	       x= (double)rand()/(double)RAND_MAX;
	       if(p >= x) {index= 1;} else {index= 0;}
	       n_pos += index;
	       i++;
	    } while (i < n);
	
		return n_pos;
	}
}

// Random number from normal distribution: It retunrs a random number according to a normal distribution ("mu", "sigma").
// CAUTION!!!!!!!!!!!!!!!!!!!! SOMETHING IS WRONG WITH THIS NORMAL DISTRIBUTION GENERATOR.
// PREV VERSION: I allowed "u" being zero, so the returned value can be -Inf.
// I fixed it, but need to test it before use it.
double _Random_normal(double mu, double sd) {
	double two_pi= 2.0 * acos(-1.0); // acos(-1) is PI
	double u, v, z;

	do {
		u= (double)rand()/(double)RAND_MAX;
		v= (double)rand()/(double)RAND_MAX;
	} while (u == 0.0);

	z= sqrt(-2.0 * log(u)) * cos(two_pi * v);

	return (mu + z * sd);
}

// Random number between two numbers: It retunrs a random number according to a uniform distribution ("a", "b").
double _Random_uniform(double a, double b) {
	double small;
	double diff= fabs(b - a);
	double x;

	if (a > b) {small= b;} else if (a < b) {small= a;} else {return a;}

	x= diff * (double)rand() / (double)(RAND_MAX/1.0);
	small= small + x;

	return small;
}

// Add a new animal to the list: It appends a node of animal to the existing list of animals.
animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list) {
	animal_node *previous_node;
	animal_node *current_node;

	previous_node= animal_list;
	current_node= animal_list;

	if (current_node == NULL) {
		animal_list= new_node;
		return (animal_list);
	} else if (current_node->demo_group == -99999) {
		new_node->next_node= current_node;
		animal_list= new_node;
		return (animal_list);
	} else {
		while (current_node->demo_group == 99999) {
			previous_node= current_node;
			current_node= current_node->next_node;
		}
		previous_node->next_node= new_node;
		new_node->next_node= current_node;
		return (animal_list);
	}
}

// Change group size based on the demographic group: It changes the number of animals in a certain group (as a combination of "demo_group" and "inf_stat").
void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals) {
	group_size[farm_id][(inf_stat * 7 + demo_group)] += n_animals;
}

// Create animal list: It generates a list of cattle for inital demographic setting of simulation.
void _Create_animal_list(int **farm_info, int **group_size, animal_node *animal_list[]) {     
	animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list);
	int _Random_binomial(int n, double p);
	double _Random_uniform(double a, double b);
    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  

	int i, j, cow_size, bh_size, yh_size;
	animal_node *new_node;

	for (i= 0; i< n_farms; i++) {
		new_node= (animal_node *)malloc(sizeof(animal_node));
		new_node->demo_group= -99999;
		new_node->dam_node= NULL;
		new_node->next_node= NULL;
		animal_list[i]= new_node;
	
		new_node= (animal_node *)malloc(sizeof(animal_node));
	    new_node->demo_group= 99999;
		new_node->dam_node= NULL;
	    new_node->next_node= NULL;
		animal_list[i]= _Add_animal_to_list(new_node, animal_list[i]);

		cow_size= farm_info[i][0];
		bh_size= (int)round((double)farm_info[i][0] * 0.22 * 0.9314);
		yh_size= (int)round((double)farm_info[i][0] * 0.22 * 0.9314); // 20% Replacement * 93.14% survival 1 - (0.2 * 0.143 + 0.8 * 0.05)
		
		// Young heifers
        for (j= 0; j< yh_size; j++) {
			new_node= (animal_node *)malloc(sizeof(animal_node));
	
			new_node->demo_group= 2;
			new_node->age= 365 - 1;
			new_node->sex= 0;
			new_node->weight= 570.0 * pow(1.0 - 0.590296 * exp(-1.0 * 0.004487 * (double)(new_node->age)), 3);
			new_node->cull_code= 0;
			new_node->days_to_cull= 99999;
			new_node->is_pregnant= 0;
			new_node->is_aied= 0;
			new_node->calf_fate= 0;
			new_node->parity= 0;
			new_node->days_to_calv= 0;
			new_node->days_to_estrous= 0;
			new_node->is_milking= 0;
			new_node->inf_stat= 0;
			new_node->days_left= 0;
			new_node->gest_at_inf= 0;
			new_node->is_vaccinated= 0;
			new_node->vac_left= 0;
			new_node->dam_node= NULL;
		    new_node->next_node= NULL;
			animal_list[i]= _Add_animal_to_list(new_node, animal_list[i]);
			_Change_group_size(group_size, i, new_node->demo_group, new_node->inf_stat, 1);  
		}

		// Breeding heifers
        for (j= 0; j< bh_size; j++) {
			new_node= (animal_node *)malloc(sizeof(animal_node));
	
			new_node->demo_group= 3;
			new_node->age= 2 * 365 - 1;
			new_node->sex= 0;
			new_node->weight= 570.0 * pow(1.0 - 0.590296 * exp(-1.0 * 0.004487 * (double)(new_node->age)), 3);
			new_node->cull_code= 0;
			new_node->days_to_cull= 99999;
			new_node->is_pregnant= 1;
			new_node->is_aied= 0;
			new_node->calf_fate= 0;
			new_node->parity= 1;
			new_node->days_to_calv= 0;
			new_node->days_to_estrous= 0;
			new_node->is_milking= 0;
			new_node->inf_stat= 0;
			new_node->days_left= 0;
			new_node->gest_at_inf= 0;
			new_node->is_vaccinated= 0;
			new_node->vac_left= 0;
			new_node->dam_node= NULL;
		    new_node->next_node= NULL;
			animal_list[i]= _Add_animal_to_list(new_node, animal_list[i]);
			_Change_group_size(group_size, i, new_node->demo_group, new_node->inf_stat, 1);  
		}

        for (j= 0; j< cow_size; j++) {
			new_node= (animal_node *)malloc(sizeof(animal_node));
	
			new_node->demo_group= 4;
			new_node->age= 3 * 365 - 1;
			new_node->sex= 0;
			new_node->weight= 570.0 * pow(1.0 - 0.590296 * exp(-1.0 * 0.004487 * (double)(new_node->age)), 3);
			new_node->cull_code= 0;
			new_node->days_to_cull= 99999;
			new_node->is_pregnant= 1;
			new_node->is_aied= 0;
			new_node->calf_fate= 0;
			new_node->parity= 2;
			new_node->days_to_calv= 0;
			new_node->days_to_estrous= 0;
			new_node->is_milking= 0;
			new_node->inf_stat= 0;
			new_node->days_left= 0;
			new_node->gest_at_inf= 0;
			new_node->is_vaccinated= 0;
			new_node->vac_left= 0;
			new_node->dam_node= NULL;
		    new_node->next_node= NULL;
			animal_list[i]= _Add_animal_to_list(new_node, animal_list[i]);
			_Change_group_size(group_size, i, new_node->demo_group, new_node->inf_stat, 1);  
		}
	}
}

// BVDV infection status of calf that will be born. SHOULD BE UPDATED WHEN DAM IS EITHER CONCEIVED OR INFECTED!!!
int _Fate_of_calf(int inf_stat, int gest_at_inf) {
	double _Random_normal(double mu, double sd);
	
	// 0: Susceptible, 1: TI, 2: Recovered, 3: iMmuned, 4: PI, 5: Vaccinated
	// From Susceptible: Susceptible calves
	if (inf_stat == 0) {
		return 0;
	// From PI: PI calves (abortion rate was assumed to be the same across TI and PI)
	} else if (inf_stat == 4) {
		return 4;
	// From Infected
	} else if (inf_stat == 1) {
		// Early pregnancy
		if (gest_at_inf >= 239) { // Day 0 ~ Day 41
			return 3; // Immuned calf from TI cattle
		// Mid pregnancy
		} else if (gest_at_inf >= 130 && gest_at_inf < 239) { // Day 42 ~ Day 150 (Pregnancy period was categorised as same as Damman et al., 2015)
			if (_Random_binomial(1, 0.933) == 1) { // 0.933 is conditional probability of calving PI (Cond: No abortion & CD)
				return 4;
			} else {
				if (_Random_binomial(1, 0.5) == 1) { // 0.50 is conditional probability of calving immuned (Cond: No abortion, No PI)
					return 3; // Immuned calf from TI cattle
				} else {
					return 2; // Recovered calf from TI cattle
				}
			}
		// Late pregnancy
		} else if (gest_at_inf >= 1 && gest_at_inf < 130) {
			return 2;
		// Infected before conception
		} else {
			return 3;
		}
	// From Recovered or immuned
	} else if (inf_stat == 2) {
		return 3;
	// From Vaccinated
	} else if (inf_stat == 5) {
	// From Infected, Recovered or immuned: Infection status of calves depends on the time of infection
		// Early pregnancy
		if (gest_at_inf >= 239) { // Day 0 ~ Day 41
			return 3; // Immuned calf from TI cattle
		// Mid pregnancy
		} else if (gest_at_inf >= 130 && gest_at_inf < 239) { // Day 42 ~ Day 150 (Pregnancy period was categorised as same as Damman et al., 2015)
			if (_Random_binomial(1, 0.2319) == 1) { // 0.2319 is conditional probability of calving PI (Cond: No abortion)
				return 4;
			} else {
				if (_Random_binomial(1, 0.7273) == 1) { // 0.7273 is conditional probability of calving immuned (Cond: No abortion, No PI)
					return 3; // Immuned calf from TI cattle
				} else {
					return 2; // Recovered calf from TI cattle
				}
			}
		// Late pregnancy
		} else if (gest_at_inf >= 1 && gest_at_inf < 130) {
			return 2;
		// Infected before conception
		} else {
			return 3;
		}
	} else {
		return 0;
	}
}

// Create a new calf: It appends a node of new calf to a list of cattle according to its dam's status.
void _Add_new_calf(int year[1], int **n_pi, int **group_size, int farm_id, animal_node *dam_node, animal_node *animal_list[]) {
	int _Random_binomial(int n, double p);
	double _Random_normal(double mu, double sd);
	animal_node *_Add_animal_to_list(animal_node *new_node, animal_node *animal_list);
    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  

	animal_node *new_node;
	
	new_node= (animal_node *)malloc(sizeof(animal_node));
	if (dam_node->demo_group == 4) {
		new_node->demo_group= 1;
	} else {
		new_node->demo_group= 0;
	}
	new_node->age= 0;
	new_node->sex= _Random_binomial(1, 0.5);
	if (new_node->sex == 0) {
		new_node->weight= _Random_normal(39.8, 0.6); // Hickson et al., 2015
	} else {
		new_node->weight= _Random_normal(43.2, 0.6);
	}
	if (dam_node->demo_group == 3) {
		new_node->weight *= 0.8811; // Hickson et al., 2012
	}
	new_node->cull_code= 0;
	new_node->days_to_cull= 99999;
	new_node->is_pregnant= 0;
	new_node->is_aied= 0;
	new_node->calf_fate= 0;
	new_node->parity= 0;
	new_node->days_to_calv= 0;
	new_node->days_to_estrous= 0;
	new_node->is_milking= 0;
	new_node->inf_stat= dam_node->calf_fate;
	new_node->gest_at_inf= 0;
	if (new_node->inf_stat == 3) {
		new_node->days_left= (int)round(_Random_normal(155.0, 31.0));
	} else {
		new_node->days_left= 0;
	}
	new_node->is_vaccinated= 0;
	new_node->vac_left= 0;
    new_node->dam_node= dam_node;
    new_node->next_node= NULL;

	if (new_node->inf_stat == 4) {
		n_pi[farm_id][year[0]] += 1;
	}

	_Change_group_size(group_size, farm_id, new_node->demo_group, new_node->inf_stat, 1);
	animal_list[farm_id]= _Add_animal_to_list(new_node, animal_list[farm_id]);
}

// Vaccination
int _Vaccination(int inf_stat, int is_vaccinated) { 
// Animal's "is_vaccinated" status should be add before running the "vaccination" function.
	double vac; 

	if (inf_stat == 0 || inf_stat == 5) { // For those susceptible or vaccinated
		if (is_vaccinated == 1) { // First vaccination: Assumed nothing happens
			return inf_stat;
//		} else if (is_vaccinated == 2) { // Second vaccination: Assumed either 1 of 3 outcomes happens
		} else if (is_vaccinated >= 2) { // Second vaccination: Assumed either 1 of 3 outcomes happens
			vac= (double)rand()/(double)RAND_MAX;
			if (vac <= 0.193) {return 0; // Not sufficient to prevent both clinical diseases and fetal infection
			} else if (vac <= 0.928) {return 5; // Sufficient to prevent clinical diseases but not fetal infection
			} else {return 2; // Sufficient to prevent both clinical diseases and fetal infection
			}
		}
//		} else {
//			return 2;
//		}
	} else { // For others than susceptible or vaccinated
		return inf_stat;
	}
}

// ELISA tests
int _ELISA(int inf_stat) {
	int _Random_binomial(int n, double p);
	double Se, Sp;
	Se= 0.835; // Lanyon et al., 2014
	Sp= 0.994; // Lanyon et al., 2014

	if (inf_stat == 1 || inf_stat == 4) {
		return (_Random_binomial(1, Se));
	} else {
		return (_Random_binomial(1, 1 - Sp));		
	}
}

// Result of test and cull
void _Test_cull(int farm_id, int **n_test, int pool, int year[1], animal_node *current_node) {
	int _Random_binomial(int n, double p);
	int _ELISA(int inf_stat);

	// Individual test
	if (pool == 0) {
		n_test[farm_id][year[0]] +=1; // First test
		if (_ELISA(current_node->inf_stat) == 1) { // Positive at the first test
			// If PI
			if (current_node->inf_stat == 4) {
				n_test[farm_id][year[0]] +=1; // Re-test after 4 weeks
				if (_ELISA(current_node->inf_stat) == 1) {
					current_node->cull_code= 1;
					current_node->days_to_cull= 28;
					if (current_node->dam_node != NULL) {
						n_test[farm_id][year[0]] +=1; // Test of the dam
						if (_ELISA(current_node->dam_node->inf_stat) == 1) {
							if (current_node->dam_node->days_to_cull > 28) {
								current_node->dam_node->cull_code= 1;
								current_node->dam_node->days_to_cull= 28;
							}
						}
					}
				}
			// If TI
			} else if (current_node->inf_stat == 1) {
				n_test[farm_id][year[0]] +=1; // Re-test after 4 weeks
				if (_ELISA(2) == 1) { // TI animals would be recovered by the second test (28 days later) 
					current_node->cull_code= 1;
					current_node->days_to_cull= 28;
					if (current_node->dam_node != NULL) {
						n_test[farm_id][year[0]] +=1; // Test of the dam
						if (_ELISA(current_node->dam_node->inf_stat) == 1) {
							if (current_node->dam_node->days_to_cull > 28) {
								current_node->dam_node->cull_code= 1;
								current_node->dam_node->days_to_cull= 28;
							}
						}
					}
				}
			// Others	
			} else {
				n_test[farm_id][year[0]] +=1; // Re-test after 4 weeks
				if (_ELISA(current_node->inf_stat) == 1) {
					current_node->cull_code= 1;
					current_node->days_to_cull= 28;
					if (current_node->dam_node != NULL) {
						n_test[farm_id][year[0]] +=1; // Test of the dam
						if (_ELISA(current_node->dam_node->inf_stat) == 1) {
							if (current_node->dam_node->days_to_cull > 28) {
								current_node->dam_node->cull_code= 1;
								current_node->dam_node->days_to_cull= 28;
							}
						}
					}
				}
			}	
		}
	// Follow-up test of previous pooled test
	} else {
		n_test[farm_id][year[0]] +=1; // Re-test after 4 weeks
		if (_ELISA(current_node->inf_stat) == 1) {
			current_node->cull_code= 1;
			current_node->days_to_cull= 0;
		}
	}
}

// Update demographic status: It updates animal's demographic compartment for every farm according to farm's demographic event date. 
void _Update_demo(int year[1], int day[1], int **n_pi, int **farm_info, int **group_size, int **n_replace, int **n_purchase,
	int **n_vaccine, int **n_test, int *ag_tested, int *follow_up, int *SCR, int *VAC, int *DF, animal_node *animal_list[]) {

    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  
	double _Random_uniform(double a, double b);
	int _Random_binomial(int n, double p);
	int _Fate_of_calf(int inf_stat, int gest_at_inf);
	void _Add_new_calf(int year[1], int **n_pi, int **group_size, int farm_id, animal_node *dam_node, animal_node *animal_list[]);
	int _Vaccination(int inf_stat, int is_vaccinated);
	int _ELISA(int inf_stat);
	void _Test_cull(int farm_id, int **n_test, int pool, int year[1], animal_node *current_node);

	int i;
	int d_psm, d_mix, d_pregscan, d_wean;
	int yr_start= yr_burn_demo + yr_burn_bvdv;
	//int yr_start= yr_burn_demo + yr_burn_bvdv - 2;
	double p_calf_mort_hfr, p_calf_mort_cow, p_conc_bull, rr_conc;
	double a_von, b_von, k_von, lwg;
	a_von= 570.0;
	b_von= 0.590296;
	k_von= 0.004487;

	p_calf_mort_hfr= 0.00086;
	p_calf_mort_cow= 0.00028;
	rr_conc= 0.6805;

	animal_node *previous_node;
	animal_node *current_node;

	for (i= 0; i< n_farms; i++) {
		d_psm= farm_info[i][1];
		d_wean= farm_info[i][2];
		d_mix= farm_info[i][3];
		d_pregscan= farm_info[i][4];

		previous_node= animal_list[i];
		current_node= animal_list[i];

		p_conc_bull= 0.61;
		lwg= 0.0;
		
		//if (day[0] == d_wean) { // Number of female calves for replacement should be reset at the day of weaning.
		if (day[0] == d_wean - 1) { // Number of female calves for replacement should be reset at the day of weaning.
			n_replace[i][0]= (int)round((double)farm_info[i][0] * 0.22);
		}
		if (day[0] == d_pregscan) {
			n_replace[i][1]= 0;
		}

		// BTM or Pooled PCR
		if ((*SCR) == 2 && year[0] >= yr_start && day[0] == 59) {ag_tested[i]= 0; follow_up[i]= 0;} // 1st March
		
		while (current_node != NULL) {
			if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
				previous_node= current_node;
				current_node= current_node->next_node;
			} else {
				// Calves_BH
				if (current_node->demo_group == 0) {
					current_node->days_to_cull -=1;
					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL; // Disconnect the pedigree line (link to dam node) before the day of dam being culled.
					current_node->age +=1;
					if (current_node->weight > 570.0) {
						lwg= 0.0;
					} else  {
						lwg= 3.0 * a_von * k_von * b_von * exp(-1.0 * k_von * (double)(current_node->age)) * pow((1.0 - b_von * exp(-1.0 * k_von * (double)(current_node->age))), 2.0);
					}
					if (current_node->inf_stat == 4) {current_node->weight += 0.810 * lwg;} else if (current_node->inf_stat == 1) {current_node->weight += 0.924 * lwg;} else {current_node->weight += lwg;}

					// Natural mortality of until weaning
					if (day[0] != d_mix) {
						if (_Random_binomial(1, p_calf_mort_hfr) == 1) {
							current_node->cull_code= 0;
							current_node->days_to_cull= 0;
						}
					// Mixing
					} else {
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						current_node->demo_group= 1;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
					}

					// Move to the next node
					previous_node= current_node;
					current_node= current_node->next_node;
				// Calves_MA
				} else if (current_node->demo_group == 1) {
					current_node->days_to_cull -=1;
					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL; // Disconnect the pedigree line (link to dam node) before the day of dam being culled.
					current_node->age +=1;
					if (current_node->weight > 570.0) {
						lwg= 0.0;
					} else  {
						lwg= 3.0 * a_von * k_von * b_von * exp(-1.0 * k_von * (double)(current_node->age)) * pow((1.0 - b_von * exp(-1.0 * k_von * (double)(current_node->age))), 2.0);
					}
					if (current_node->inf_stat == 4) {current_node->weight += 0.810 * lwg;} else if (current_node->inf_stat == 1) {current_node->weight += 0.924 * lwg;} else {current_node->weight += lwg;}

					// Natural mortality of until weaning
					if (day[0] != d_wean) {
						if (_Random_binomial(1, p_calf_mort_cow) == 1) {
							current_node->cull_code= 0;
							current_node->days_to_cull= 0;
						}
					// Weaning
					} else {
						// Selection as YHs
						if (current_node->sex == 0 && current_node->age >= 173) {
							if (n_replace[i][0] > 0) {
								_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
								current_node->demo_group= 2;

								// Test and cull
								if ((*SCR) == 1 && year[0] >= yr_start) { // 1st March
									_Test_cull(i, n_test, 0, year, current_node);
								}

								// Follow-up test
								if ((*SCR) == 2 && year[0] >= yr_start && day[0] == d_wean && follow_up[i] > 0) {
									_Test_cull(i, n_test, 1, year, current_node);
								}
					
								// Vaccination: 1st vaccination for VP 1
								if (((*VAC) == 1 || (*VAC) == 2) && year[0] >= yr_start) {
									if (day[0] == d_wean) { // First vaccination
										n_vaccine[i][year[0]] +=1;
										current_node->is_vaccinated= 1;
										current_node->inf_stat= _Vaccination(current_node->inf_stat, current_node->is_vaccinated);
									}
								}

								_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
								n_replace[i][0] -=1;
							// Surplus as heifers/stores
							} else {
								_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
								current_node->demo_group= 5;
								_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
							}
						} else {
							// Surplus as bulls/steers
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->demo_group= 6;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						}
					}

					// Move to the next node
					previous_node= current_node;
					current_node= current_node->next_node;
				// Young heifers
				} else if (current_node->demo_group == 2) {
					current_node->days_to_cull -=1;
					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL; // Disconnect the pedigree line (link to dam node) before the day of dam being culled.
					current_node->age +=1;
					if (current_node->weight > 570.0) {
						lwg= 0.0;
					} else  {
						lwg= 3.0 * a_von * k_von * b_von * exp(-1.0 * k_von * (double)(current_node->age)) * pow((1.0 - b_von * exp(-1.0 * k_von * (double)(current_node->age))), 2.0);
					}
					if (current_node->inf_stat == 4) {current_node->weight += 0.810 * lwg;} else if (current_node->inf_stat == 1) {current_node->weight += 0.924 * lwg;} else {current_node->weight += lwg;}

					// Screening
					if (((*SCR) == 2) && day[0] == 59 && year[0] >= yr_start) { // 1st March
						if (ag_tested[i] < 15) {
							ag_tested[i] +=1;
							n_test[i][year[0]] +=1;
							if (_ELISA(current_node->inf_stat) == 1) {
								follow_up[i] +=1;
							}
						}
					}
					// Follow-up test
					if ((*SCR) == 2 && year[0] >= yr_start && day[0] == d_wean && follow_up[i] > 0) {
						_Test_cull(i, n_test, 1, year, current_node);
					}
					
					// Vaccination: 2nd & 3rd vaccination for VP 1
					if (((*VAC) == 1 || (*VAC) == 2) && year[0] >= yr_start) {
						// 2nd vaccination: May 1st
						if (day[0] == 120) {
							n_vaccine[i][year[0]] +=1;
							current_node->is_vaccinated +=1;
							current_node->vac_left= 180;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= _Vaccination(current_node->inf_stat, current_node->is_vaccinated);
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						}
					}
					// 3rd vaccination: October 29th
					if (((*VAC) == 1 || (*VAC) == 2) && year[0] > yr_start) {
						if (day[0] == 301) {
							n_vaccine[i][year[0]] +=1;
							current_node->is_vaccinated +=1;
							current_node->vac_left= 365;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= _Vaccination(current_node->inf_stat, current_node->is_vaccinated);
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						}
					} else if (((*VAC) == 1 || (*VAC) == 2) && year[0] == yr_start) {
						if (day[0] == 301) {
							n_vaccine[i][year[0]] +=3;
							current_node->is_vaccinated +=3;
							current_node->vac_left= 365;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= _Vaccination(current_node->inf_stat, current_node->is_vaccinated);
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						}
					}
					// Vaccination: 1st vaccination for VP 2
					if (((*VAC) == 3 || (*VAC) == 4) && year[0] >= yr_start) {
						if (day[0] == 273) { // First vaccination
							n_vaccine[i][year[0]] +=1;
							current_node->is_vaccinated= 1;
							current_node->inf_stat= _Vaccination(current_node->inf_stat, current_node->is_vaccinated);
						}
					}

					// At the puberty
					if (current_node->age == 380) { // Puberty starts at 380 days old. McNaughton et al., 2002
						current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // Set the counter upto the second heat (silent first heat)
					// After the puberty
					} else if (current_node->age > 380) {
						// Pregnant heifers
						if (current_node->is_pregnant == 1) {
							// Natural abortion
							if (_Random_binomial(1, 0.00013) == 1) {
								current_node->is_pregnant= 0;
								current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0));// Set the counter upto the second heat (womb preparation + silent first heat)
								current_node->days_to_calv= 0;
								current_node->gest_at_inf= 0;
							} else {
								// Calving
								if (current_node->days_to_calv == 0) {
									_Add_new_calf(year, n_pi, group_size, i, current_node, animal_list);
									current_node->parity +=1;
									current_node->is_pregnant= 0;
									current_node->is_milking= 1;
									current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0)); 
									// Set the counter upto the second heat (womb preparation (calving to first oestrous (Olds et al., 1953) + silent first heat)
									current_node->gest_at_inf= 0;
								} else {
									current_node->days_to_calv -=1; // REDUCE THE DAYS TO CALV
								}
							}
						// Empty heifers
						} else {
							// On Heat
							if (current_node->days_to_estrous == 0) {
								// DURING THE 9 WEEKS OF BREEDING PERIOD
								if ((365 + day[0] - d_psm) % 365 >= 0 && (365 + day[0] - d_psm) % 365 < 7 * 9) {
									// TI or PI heifers
									if (current_node->inf_stat == 1 || current_node->inf_stat == 4) {
										if (_Random_binomial(1, p_conc_bull * rr_conc) == 1) {
											current_node->is_pregnant= 1;
											current_node->calf_fate= _Fate_of_calf(current_node->inf_stat, current_node->gest_at_inf);
											current_node->days_to_calv= 280; // It means the calving will occurr after 281 days as calving occurs under the condition of days_to_calv = 0.
										} else {
											current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
										}
									} else {
										if (_Random_binomial(1, p_conc_bull) == 1) {
											//printf("Conc.\n");
											current_node->is_pregnant= 1;
											current_node->calf_fate= _Fate_of_calf(current_node->inf_stat, current_node->gest_at_inf);
											current_node->days_to_calv= 280; // It means the calving will occurr after 281 days as calving occurs under the condition of days_to_calv = 0.
										} else {
											current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
										}
									}
								// Out of mating period
								} else {
									current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
								}
							// No heat
							} else {
								current_node->days_to_estrous -=1; // REDUCE THE DAYS TO ESTROUS
							}
						}
					}

					// Becoming BHs
					if (day[0] == d_wean) {
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						current_node->demo_group= 3;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
					}
					
					// Move to the next node
					previous_node= current_node;
					current_node= current_node->next_node;
				// Breeding heifers
				} else if (current_node->demo_group == 3) {
					current_node->days_to_cull -=1;
					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL; // Disconnect the pedigree line (link to dam node) before the day of dam being culled.
					current_node->age +=1;
					if (current_node->weight > 570.0) {
						lwg= 0.0;
					} else  {
						lwg= 3.0 * a_von * k_von * b_von * exp(-1.0 * k_von * (double)(current_node->age)) * pow((1.0 - b_von * exp(-1.0 * k_von * (double)(current_node->age))), 2.0);
					}
					if (current_node->inf_stat == 4) {current_node->weight += 0.810 * lwg;} else if (current_node->inf_stat == 1) {current_node->weight += 0.924 * lwg;} else {current_node->weight += lwg;}

					// Follow-up test
					if ((*SCR) == 2 && year[0] >= yr_start && day[0] == d_wean && follow_up[i] > 0) {
						_Test_cull(i, n_test, 1, year, current_node);
					}

					// Vaccination: Annual booster for VP 1
					if ((*VAC) == 2 && year[0] > yr_start) {
						if (day[0] == 301) { // October 29th
							n_vaccine[i][year[0]] +=1;
							current_node->is_vaccinated +=1;
							current_node->vac_left= 365;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= _Vaccination(current_node->inf_stat, current_node->is_vaccinated);
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						}
					} else if ((*VAC) == 2 && year[0] == yr_start) {
						if (day[0] == 301) { // October 29th
							n_vaccine[i][year[0]] +=3;
							current_node->is_vaccinated +=3;
							current_node->vac_left= 365;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= _Vaccination(current_node->inf_stat, current_node->is_vaccinated);
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						}
					}
					// Vaccination: 2nd vaccination for VP 2
					if (((*VAC) == 3 || (*VAC) == 4) && year[0] >= yr_start) {
						if (day[0] == 301) {
							n_vaccine[i][year[0]] +=1;
							current_node->is_vaccinated +=1;
							current_node->vac_left= 180;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= _Vaccination(current_node->inf_stat, current_node->is_vaccinated);
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						}
					}
				
					// Pregnant heifers
					if (current_node->is_pregnant == 1) {
						// Natural abortion
						if (_Random_binomial(1, 0.00013) == 1) {
							current_node->is_pregnant= 0;
							current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0));// Set the counter upto the second heat (womb preparation + silent first heat)
							current_node->days_to_calv= 0;
							current_node->gest_at_inf= 0;
						} else {
							// Calving
							if (current_node->days_to_calv == 0) {
								_Add_new_calf(year, n_pi, group_size, i, current_node, animal_list);
								current_node->parity +=1;
								current_node->is_pregnant= 0;
								current_node->is_milking= 1;
								current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0)); 
								// Set the counter upto the second heat (womb preparation (calving to first oestrous (Olds et al., 1953) + silent first heat)
								current_node->gest_at_inf= 0;
							} else {
								current_node->days_to_calv -=1; // REDUCE THE DAYS TO CALV
							}
						}
					// Empty heifers
					} else {
						// On Heat
						if (current_node->days_to_estrous == 0) {
							// DURING THE 9 WEEKS OF BREEDING PERIOD
							if ((365 + day[0] - d_psm) % 365 >= 0 && (365 + day[0] - d_psm) % 365 < 7 * 9) {
								// TI or PI heifers
								if (current_node->inf_stat == 1 || current_node->inf_stat == 4) {
									if (_Random_binomial(1, p_conc_bull * rr_conc) == 1) {
										current_node->is_pregnant= 1;
										current_node->calf_fate= _Fate_of_calf(current_node->inf_stat, current_node->gest_at_inf);
										current_node->days_to_calv= 280; // It means the calving will occurr after 281 days as calving occurs under the condition of days_to_calv = 0.
									} else {
										current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
									}
								} else {
									if (_Random_binomial(1, p_conc_bull) == 1) {
										current_node->is_pregnant= 1;
										current_node->calf_fate= _Fate_of_calf(current_node->inf_stat, current_node->gest_at_inf);
										current_node->days_to_calv= 280; // It means the calving will occurr after 281 days as calving occurs under the condition of days_to_calv = 0.
									} else {
										current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
									}
								}
							// Out of mating period
							} else {
								current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
							}
						// No heat
						} else {
							current_node->days_to_estrous -=1; // REDUCE THE DAYS TO ESTROUS
						}
					}

					// Becoming MAs
					if (day[0] == d_mix) {
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						current_node->demo_group= 4;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
					}

					// Move to the next node
					previous_node= current_node;
					current_node= current_node->next_node;
				// MA cows
				} else if (current_node->demo_group == 4) {
					current_node->days_to_cull -=1;
					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL; // Disconnect the pedigree line (link to dam node) before the day of dam being culled.
					current_node->age +=1;
					if (current_node->weight > 570.0) {
						lwg= 0.0;
					} else  {
						lwg= 3.0 * a_von * k_von * b_von * exp(-1.0 * k_von * (double)(current_node->age)) * pow((1.0 - b_von * exp(-1.0 * k_von * (double)(current_node->age))), 2.0);
					}
					if (current_node->inf_stat == 4) {current_node->weight += 0.810 * lwg;} else if (current_node->inf_stat == 1) {current_node->weight += 0.924 * lwg;} else {current_node->weight += lwg;}

					// Follow-up test
					if ((*SCR) == 2 && year[0] >= yr_start && day[0] == d_wean && follow_up[i] > 0) {
						_Test_cull(i, n_test, 1, year, current_node);
					}

					// Vaccination: Annual booster for VP 1 & VP 2
					if (((*VAC) == 2 || (*VAC) == 4) && year[0] > yr_start) {
						if (day[0] == 301) {
							n_vaccine[i][year[0]] +=1;
							current_node->is_vaccinated +=1;
							current_node->vac_left= 365;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= _Vaccination(current_node->inf_stat, current_node->is_vaccinated);
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						}
					} else if (((*VAC) == 2 || (*VAC) == 4) && year[0] == yr_start) {
						if (day[0] == 301) {
							n_vaccine[i][year[0]] +=3;
							current_node->is_vaccinated +=3;
							current_node->vac_left= 365;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= _Vaccination(current_node->inf_stat, current_node->is_vaccinated);
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						}
					}

					// Pregnant cows
					if (current_node->is_pregnant == 1) {
						// Natural abortion
						if (_Random_binomial(1, 0.00013) == 1) {
							current_node->is_pregnant= 0;
							current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0));// Set the counter upto the second heat (womb preparation + silent first heat)
							current_node->days_to_calv= 0;
							current_node->gest_at_inf= 0;
						} else {
							// Calving
							if (current_node->days_to_calv == 0) {
								_Add_new_calf(year, n_pi, group_size, i, current_node, animal_list);
								current_node->parity +=1;
								current_node->is_pregnant= 0;
								current_node->is_milking= 1;
								current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0)); 
								// Set the counter upto the second heat (womb preparation (calving to first oestrous (Olds et al., 1953) + silent first heat)
								current_node->gest_at_inf= 0;
							} else {
								current_node->days_to_calv -=1; // REDUCE THE DAYS TO CALV
							}
						}
					// Empty cows
					} else {
						// On Heat
						if (current_node->days_to_estrous == 0) {
							// DURING THE 9 WEEKS OF BREEDING PERIOD
							if ((365 + day[0] - d_psm) % 365 >= 0 && (365 + day[0] - d_psm) % 365 < 7 * 9) {
								// TI or PI heifers
								if (current_node->inf_stat == 1 || current_node->inf_stat == 4) {
									if (_Random_binomial(1, p_conc_bull * rr_conc) == 1) {
										current_node->is_pregnant= 1;
										current_node->calf_fate= _Fate_of_calf(current_node->inf_stat, current_node->gest_at_inf);
										current_node->days_to_calv= 280; // It means the calving will occurr after 281 days as calving occurs under the condition of days_to_calv = 0.
									} else {
										current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
									}
								} else {
									if (_Random_binomial(1, p_conc_bull) == 1) {
										current_node->is_pregnant= 1;
										current_node->calf_fate= _Fate_of_calf(current_node->inf_stat, current_node->gest_at_inf);
										current_node->days_to_calv= 280; // It means the calving will occurr after 281 days as calving occurs under the condition of days_to_calv = 0.
									} else {
										current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
									}
								}
							// Out of mating period
							} else {
								current_node->days_to_estrous= (int)round(_Random_uniform(18.0, 24.0)); // RESET THE DAYS TO ESTROUS DEPENDING ON THE PREGNANT STATUS
							}
						// No heat
						} else {
							current_node->days_to_estrous -=1; // REDUCE THE DAYS TO ESTROUS
						}
					}

					// Pregnancy scan
					if (day[0] == d_pregscan) {
						if (current_node->is_pregnant == 1) {
							n_replace[i][1] +=1;
						} else {
							// Tagged empty MA cows
							current_node->cull_code= 2;
							current_node->days_to_cull= (365 + d_wean - d_pregscan) % 365;
						}
					}
	
					// Fate of MAs
					if (day[0] == d_wean) {
						// In case of the pregnant MAs > target MA cow herd size: MAs are culled whether it's pregnant or not.
						if (n_replace[i][1] > farm_info[i][0]) {
							// Pregnant MAs are removed.
							if (year[0] < yr_burn_demo) {
								if (current_node->cull_code == 0 && current_node->parity >= 2) {
									current_node->cull_code= 3;
									current_node->days_to_cull= 0;
									n_replace[i][1] -=1;						
								}
							} else {
								if (current_node->cull_code == 0 && current_node->parity >= 6) {
									current_node->cull_code= 3;
									current_node->days_to_cull= 0;
									n_replace[i][1] -=1;						
								}
							}
						// In case of the pregnant MAs < target MA cow herd size: Some of empty MAs survive.
						} else if (n_replace[i][1] < farm_info[i][0]) {
							// Empty MAs are redempted.
							if (current_node->cull_code == 2) {
								current_node->cull_code= 0;
								current_node->days_to_cull= 99999;
								current_node->parity= 0;
								current_node->inf_stat= 0;
								//if (_Random_binomial(1, 0.02) == 1) {
								//	current_node->inf_stat= 4;
								//} else {
								//	current_node->inf_stat= 0;
								//}
								current_node->is_vaccinated= 0;
								n_purchase[i][year[0]] +=1;
								n_replace[i][1] +=1;
							}
						}
						// Do nothing when the pregnant MAs = target MAs
					}

					// Move to the next node
					previous_node= current_node;
					current_node= current_node->next_node;
				// Heifers/Stores
				} else if (current_node->demo_group == 5) {
					current_node->days_to_cull -=1;
					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL; // Disconnect the pedigree line (link to dam node) before the day of dam being culled.
					current_node->age +=1;
					if (current_node->weight > 570.0) {
						lwg= 0.0;
					} else  {
						lwg= 3.0 * a_von * k_von * b_von * exp(-1.0 * k_von * (double)(current_node->age)) * pow((1.0 - b_von * exp(-1.0 * k_von * (double)(current_node->age))), 2.0);
					}
					if (current_node->inf_stat == 4) {current_node->weight += 0.810 * lwg;} else if (current_node->inf_stat == 1) {current_node->weight += 0.924 * lwg;} else {current_node->weight += lwg;}

					// TO ABATTOIR
					if (current_node->age >= 820) {
						current_node->cull_code= 2;
						current_node->days_to_cull= 0;
					}

					// Move to the next node
					previous_node= current_node;
					current_node= current_node->next_node;
				// Bulls/Steers	
				} else if (current_node->demo_group == 6) {
					current_node->days_to_cull -=1;
					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL; // Disconnect the pedigree line (link to dam node) before the day of dam being culled.
					current_node->age +=1;
					if (current_node->weight > 570.0) {
						lwg= 0.0;
					} else  {
						lwg= 3.0 * a_von * k_von * b_von * exp(-1.0 * k_von * (double)(current_node->age)) * pow((1.0 - b_von * exp(-1.0 * k_von * (double)(current_node->age))), 2.0);
					}
					if (current_node->inf_stat == 4) {current_node->weight += 0.810 * lwg;} else if (current_node->inf_stat == 1) {current_node->weight += 0.924 * lwg;} else {current_node->weight += lwg;}

					// TO ABATTOIR
					if (current_node->age >= 670) {
						current_node->cull_code= 2;
						current_node->days_to_cull= 0;
					}

					// Move to the next node
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			}
		}
	}
}

// Cull/sell cattle based on the culling status
void _Update_cull(int year[1], int **group_size, double **wgt_beef, int **n_surplus, animal_node *animal_list[]) {
	int i;

	animal_node *previous_node;
	animal_node *current_node;

	for (i= 0; i< n_farms; i++) {
		previous_node= animal_list[i];
		current_node= animal_list[i];

		while (current_node != NULL) {
			if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
				previous_node= current_node;
				current_node= current_node->next_node;
			} else {
				// Calves_BH or _MA
				if (current_node->demo_group == 0 || current_node->demo_group == 1) {
					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL;
					if (current_node->days_to_cull == 0) {
						if (current_node->cull_code != 0) {wgt_beef[i][year[0]] += 0.5 * 0.5990 * current_node->weight;} // test pos. calves
						previous_node->next_node= current_node->next_node;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						free(current_node);
						current_node= previous_node->next_node;
					} else {
						previous_node= current_node;
						current_node= current_node->next_node;
					}
				// Young heifers
				} else if (current_node->demo_group == 2) {
					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL;
					if (current_node->days_to_cull == 0) {
						if (current_node->cull_code != 0) {wgt_beef[i][year[0]] += 0.5 * 0.5990 * current_node->weight;} // test pos. YHs
						previous_node->next_node= current_node->next_node;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						free(current_node);
						current_node= previous_node->next_node;
					} else {
						previous_node= current_node;
						current_node= current_node->next_node;
					}
				// Breeding heifers
				} else if (current_node->demo_group == 3) {
 					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL;
					if (current_node->days_to_cull == 0) {
						if (current_node->cull_code != 0) { // empty or test pos. BH
							wgt_beef[i][year[0]] += 0.5 * current_node->weight;
						}
						previous_node->next_node= current_node->next_node;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						free(current_node);
						current_node= previous_node->next_node;
					} else {
						previous_node= current_node;
						current_node= current_node->next_node;
					}
				// MAs
				} else if (current_node->demo_group == 4) {
					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL;
					if (current_node->days_to_cull == 0) {
						if (current_node->cull_code == 1 || current_node->cull_code == 2) { // // either empty cows or test pos. cows
							wgt_beef[i][year[0]] += 0.5 * 0.7656 * current_node->weight; 
						} else { // pregnant cows for sale
							n_surplus[i][year[0]] +=1;
						} 
						previous_node->next_node= current_node->next_node;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						free(current_node);
						current_node= previous_node->next_node;
					} else {
						previous_node= current_node;
						current_node= current_node->next_node;
					}
				// Heifers/Stores
				} else if (current_node->demo_group == 5) {
 					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL;
					if (current_node->days_to_cull == 0) {
						if (current_node->cull_code != 0) {wgt_beef[i][year[0]] += 0.5 * current_node->weight;}
						previous_node->next_node= current_node->next_node;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						free(current_node);
						current_node= previous_node->next_node;
					} else {
						previous_node= current_node;
						current_node= current_node->next_node;
					}
				// Steers/Bulls
				} else if (current_node->demo_group == 6) {
 					if (current_node->dam_node != NULL && current_node->dam_node->days_to_cull == 0) current_node->dam_node= NULL;
					if (current_node->days_to_cull == 0) {
						if (current_node->cull_code != 0) {wgt_beef[i][year[0]] += 0.5 * current_node->weight;}
						previous_node->next_node= current_node->next_node;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						free(current_node);
						current_node= previous_node->next_node;
					} else {
						previous_node= current_node;
						current_node= current_node->next_node;
					}
				}
			}
		}
	}
}

void _Chage_seropositive(double *prop, int **group_size, animal_node *animal_list[]) {
    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  

	int i;

	animal_node *previous_node;
	animal_node *current_node;

	for (i= 0; i< n_farms; i++) {
		previous_node= animal_list[i];
		current_node= animal_list[i];

		while (current_node != NULL) {
			if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
				previous_node= current_node;
				current_node= current_node->next_node;
			} else {
				if (_Random_binomial(1, (*prop)) == 1) {
					_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
					current_node->inf_stat= 2;
					_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
				}

				previous_node= current_node;
				current_node= current_node->next_node;
			}
		}

/*
		r_obs= group_size[i][18];
		r_exp= (int)round((*prop) * (double)n_ma);
		diff= r_exp - r_obs; //printf("%i.\n", diff);

		if (diff > 0) {
			while (current_node != NULL && diff > 0) {
				if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
					previous_node= current_node;
					current_node= current_node->next_node;
				} else if (current_node->inf_stat == 2 && current_node->demo_group == 4) {
					_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
					current_node->inf_stat= 0;
					_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
					diff -=1;

					previous_node= current_node;
					current_node= current_node->next_node;
				} else {
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			}
		} else if (diff < 0) {
			while (current_node != NULL && diff < 0) {
				if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
					previous_node= current_node;
					current_node= current_node->next_node;
				} else if (current_node->inf_stat == 0 && current_node->demo_group == 4) {
					_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
					current_node->inf_stat= 2;
					_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
					diff +=1;

					previous_node= current_node;
					current_node= current_node->next_node;
				} else {
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			}
		}
*/
	}
}

// BVDV introduction via purchase: Randomly designating X new born calves as PI (= Trojan cows).
void _Seed_infection(int year[1], int **group_size, int *n_troj, animal_node *animal_list[]) {
    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  

	int i;

	animal_node *previous_node;
	animal_node *current_node;

	for (i= 0; i< n_farms; i++) {
		previous_node= animal_list[i];
		current_node= animal_list[i];

		while (current_node != NULL) {
			if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
				previous_node= current_node;
				current_node= current_node->next_node;
			} else if ((current_node->demo_group == 0 || current_node->demo_group == 1) && current_node->age == 0 && current_node->sex == 0) {
				if (current_node->inf_stat != 4) {
					if (n_troj[i] > 0) {
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						current_node->inf_stat= 4;
						n_troj[i] -=1;
						current_node->days_left= 0;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);

						// Dam of the calf
						_Change_group_size(group_size, i, current_node->dam_node->demo_group, current_node->dam_node->inf_stat, -1);
						if (current_node->dam_node != NULL) {
							if (current_node->dam_node->inf_stat == 0) {
								current_node->dam_node->inf_stat= 2;
							} else if (current_node->dam_node->inf_stat == 5) {
								current_node->dam_node->inf_stat= 2;
								current_node->dam_node->is_vaccinated= 0;
							}
						}
						_Change_group_size(group_size, i, current_node->dam_node->demo_group, current_node->dam_node->inf_stat, 1);

						previous_node= current_node;
						current_node= current_node->next_node;
					} else {
						previous_node= current_node;
						current_node= current_node->next_node;
					}
				} else {
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			} else {
				previous_node= current_node;
				current_node= current_node->next_node;
			}
		}
	}
}

// Abortion: It returns an abortion indicater according to the gestation period (= "days_to_calv") of infected cattle. We assumed that there was no difference of abortion rate between TI and PI.
int _Abortion(int inf_stat, int days_to_calv) {
	int _Random_binomial(int n, double p);

	// TI
	if (inf_stat == 1) {
		// Early pregnancy (abortion probability= 0.1212)
		if (days_to_calv >= 239) {
			if (_Random_binomial(1, 0.1212) == 1) {return 1;
			} else {return 0;}
		// Mid pregnancy (abortion probability= 0.1194, congenital deformality= 0.0625)
		} else if (days_to_calv >= 130 && days_to_calv < 239) {
			if (_Random_binomial(1, 0.1744) == 1) {return 1;
			} else {return 0;}
		// Late pregnancy: No abortion
		} else {return 0;}
	// Recovered from BVDV infection due to insufficient vaccination for fetal infection
	} else if (inf_stat == 2) { 
		// Early pregnancy (abortion probability= 0.0773)
		if (days_to_calv >= 239) {
			if (_Random_binomial(1, 0.0773) == 1) {return 1;
			} else {return 0;}
		// Mid pregnancy (abortion probability= 0.0762, congenital deformality: ignored)
		} else if (days_to_calv >= 130 && days_to_calv < 239) {
			if (_Random_binomial(1, 0.0762) == 1) {return 1;
			} else {return 0;}
		// Late pregnancy: No abortion
		} else {return 0;}
	// PI
	} else if (inf_stat == 4) {
		// Early pregnancy (abortion probability= 0.1212: abortion probability/day= 0.0086)
		if (days_to_calv >= 239) {
			if (_Random_binomial(1, 0.0086) == 1) {return 1;
			} else {return 0;}
		// Mid pregnancy (abortion probability= 0.1194, congenital deformality= 0.0625: abortion+CD probability/day= 0.0128)
		} else if (days_to_calv >= 130 && days_to_calv < 239) {
			if (_Random_binomial(1, 0.0128) == 1) {return 1;
			} else {return 0;}
		// Late pregnancy: No abortion
		} else {return 0;}
	} else {return 0;
	}
}

// Disease spread: It updates BVD infection status of each cattle.
void _Update_infection(int year[1], int **farm_info, int **group_size, int **n_replace, double *betas, double *p_pi_neighbour, double *lambda, int **n_ti, int *DF, animal_node *animal_list[]) {
	// 0: Calves (F and M), 1: RHs1 (weaning to off-grazing), 2: RHs2 (separately managed), 3: MAs, 4: R1 bulls, 5: R2 bulls, 6: R3 bulls.
	int _Random_binomial(int n, double p);
	double _Random_uniform(double a, double b);
    void _Change_group_size(int **group_size, int farm_id, int demo_group, int inf_stat, int n_animals);  
	int _Abortion(int inf_stat, int days_to_calv);
	
	// Variables for simulation
	int i, j, d_psm, n_demo_group= 7;
	int n_group[n_demo_group];
	double total_prev_pi;
	double prev_pi[n_demo_group], prev_ti[n_demo_group], density_pi[n_demo_group], p_inf[n_demo_group];

	double beta_pi= betas[0]; // Beta for within-herd PI
	double beta_ti= beta_pi * betas[1]; // %Beta for within-herd TI
	double beta_bh= beta_pi * betas[2]; // %Beta for between-herd PI
	double beta_bf= beta_pi * betas[3]; // %Beta for between-farm PI
	double lambda_pi= (*lambda);
	// Considering the life-expectancy of cattle, recovered cattle were assumed to be life-long immune to BVDV infection.

	if ((*DF) == 1) {
		beta_bf= 0.0;
	}

	animal_node *previous_node;
	animal_node *current_node;

	for (i= 0; i< n_farms; i++) {
		d_psm= farm_info[i][1];
		total_prev_pi= 0.0;

		// Group sizes
		for (j= 0; j< n_demo_group; j++) {
			if (j == 0 || j == 3) {
				n_group[j]= group_size[i][0] + group_size[i][7] + group_size[i][14] + group_size[i][21] + group_size[i][28] + group_size[i][35] +
							group_size[i][3] + group_size[i][10] + group_size[i][17] + group_size[i][24] + group_size[i][31] + group_size[i][38];
			} else if (j == 1 || j == 4) {
				n_group[j]= group_size[i][1] + group_size[i][8] + group_size[i][15] + group_size[i][22] + group_size[i][29] + group_size[i][36] +
							group_size[i][4] + group_size[i][11] + group_size[i][18] + group_size[i][25] + group_size[i][32] + group_size[i][39];
			} else {
				n_group[j]= group_size[i][j] + group_size[i][j+7] + group_size[i][j+14] + group_size[i][j+21] + group_size[i][j+28] + group_size[i][j+35];
			}
		}

		// Add the number of bulls to the group size of R2 heifers & MA cows
//		if ((365 + day[0] - d_psm) % 365 >= 0 && (365 + day[0] - d_psm) % 365 < 7 * 9) {
//			n_group[3] += (int)(round((double)n_group[3] / 30)); // Hickson et al., 2012 (1:30)
//			n_group[4] += (int)(round((double)n_group[4] / 30));
//		}
	
		// Prevalence of PI and TI
		for (j= 0; j< n_demo_group; j++) {
			if (j == 0 || j == 3) {
				if (n_group[j] == 0) {
					prev_pi[j]= 0.0; prev_ti[j]= 0.0;
				} else {
					prev_pi[j]= (double)(group_size[i][28] + group_size[i][31])/n_group[j]; prev_ti[j]= (double)(group_size[i][7] + group_size[i][10])/n_group[j];
				}
				total_prev_pi += prev_pi[j] / 2.0;
			} else if (j == 1 || j == 4) {
				if (n_group[j] == 0) {
					prev_pi[j]= 0.0; prev_ti[j]= 0.0;
				} else {
					prev_pi[j]= (double)(group_size[i][29] + group_size[i][32])/n_group[j]; prev_ti[j]= (double)(group_size[i][8] + group_size[i][11])/n_group[j];
				}
				total_prev_pi += prev_pi[j] / 2.0;
			} else {
				if (n_group[j] == 0) {prev_pi[j]= 0.0; prev_ti[j]= 0.0;} else {prev_pi[j]= (double)group_size[i][j+28]/n_group[j]; prev_ti[j]= (double)group_size[i][j+7]/n_group[j];}
				total_prev_pi += prev_pi[j];
			}
		}
		
		// Between-herd infectious pressure
		for (j= 0; j< n_demo_group; j++) {
			//if (n_group[j] == 0) {density_pi[j]= 0.0;} else {density_pi[j]= (total_prev_pi - prev_pi[j])/n_group[j];}
			if (n_group[j] == 0) {density_pi[j]= 0.0;} else {density_pi[j]= (total_prev_pi - prev_pi[j])/(double)n_group[j];}
		}

		// Probabilities of infection per demographic group
		for (j= 0; j< n_demo_group; j++) {
			p_inf[j]= 1 - exp(-1.0 * (beta_pi * prev_pi[j] + beta_bh * density_pi[j] + beta_ti * prev_ti[j] + beta_bf * (*p_pi_neighbour)/(double)n_group[j]));
		}

		previous_node= animal_list[i];
		current_node= animal_list[i];

		// Update transmission
		while (current_node != NULL) {
			// 0: Susceptible, 1: TI, 2: Recovered, 3: iMmuned, 4: PI, 5: Vaccinated.
			// Recovered cattle
			if (current_node->demo_group == 99999 || current_node->demo_group == -99999) {
				previous_node= current_node;
				current_node= current_node->next_node;
			} else {
				// Susceptible -> TI
				if (current_node->inf_stat == 0) {
					if (_Random_binomial(1, p_inf[current_node->demo_group]) == 1) {
						n_ti[i][year[0]] +=1;
						// Update infection status
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						current_node->inf_stat= 1;
						current_node->is_vaccinated= 0;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
						current_node->days_left= (int)round(_Random_uniform(10.0, 20.0));
						// Reproduction outcome
						if (current_node->is_pregnant == 1) {
							if (_Abortion(current_node->inf_stat, current_node->days_to_calv) == 1) {
								current_node->is_pregnant= 0;
								current_node->days_to_calv= 0;
								current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0)); 
								// Set the counter upto the second heat (womb preparation (calving to first oestrous (Olds et al., 1953) + silent first heat)
								current_node->gest_at_inf= 0;
							} else {
								current_node->gest_at_inf= current_node->days_to_calv;
								current_node->calf_fate= _Fate_of_calf(current_node->inf_stat, current_node->gest_at_inf);
							}
						}
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				// TI -> Recovered
				} else if (current_node->inf_stat == 1) {
					// Recovery from TI
					if (current_node->days_left == 0) {
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						current_node->inf_stat= 2;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
					// Not yet recovered
					} else {
						current_node->days_left -=1;
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				// Vaccinated (full)
				} else if (current_node->inf_stat == 2) {
					if (current_node->is_vaccinated > 0) {
						current_node->vac_left -=1;
						if (current_node->vac_left <= 0) {
							current_node->inf_stat= 0;
						}
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				// Waning-off of Maternal Ab
				} else if (current_node->inf_stat == 3) {
					// iMmuned -> Susceptible
					if (current_node->days_left == 0) {
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
						current_node->inf_stat= 0;
						_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
					} else if (current_node->days_left > 0 && current_node->age >= 14) {
						// Anamnestic response
						if (_Random_binomial(1, p_inf[current_node->demo_group]) == 1) {
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= 2;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
							current_node->days_left= 0;
						} else {
							// Depletion of MAb
							current_node->days_left -=1;
						}
					} else {
						// Depletion of MAb
						current_node->days_left -=1;
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				// Mortality of PIs
				} else if (current_node->inf_stat == 4) {
					// Death of PI
					if (_Random_binomial(1, 1 - lambda_pi) == 0) {
						current_node->cull_code= 0;
						current_node->days_to_cull= 0;
					// PI survived
					} else {
						// Reproduction outcome
						if (current_node->is_pregnant == 1 && _Abortion(current_node->inf_stat, current_node->days_to_calv) == 1) {
							current_node->is_pregnant= 0;
							current_node->days_to_calv= 0;
							current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0));
							// Set the counter upto the second heat (womb preparation (calving to first oestrous (Olds et al., 1953) + silent first heat)
							current_node->gest_at_inf= 0;
						}
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				// Vaccinated
				} else if (current_node->inf_stat == 5) {
					current_node->vac_left -=1;
					if (current_node->vac_left == 0) {
						current_node->inf_stat= 0;
					} else {
						// Anamnestic response
						if (_Random_binomial(1, p_inf[current_node->demo_group]) == 1) {
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, -1);
							current_node->inf_stat= 2;
							current_node->is_vaccinated= 0;
							_Change_group_size(group_size, i, current_node->demo_group, current_node->inf_stat, 1);
							// Reproduction outcome
							if (current_node->is_pregnant == 1) {
								if (_Abortion(current_node->inf_stat, current_node->days_to_calv) == 1) {
									current_node->is_pregnant= 0;
									current_node->days_to_calv= 0;
									current_node->days_to_estrous= (int)round(_Random_uniform(15.0, 49.0) + _Random_uniform(18.0, 24.0)); 
									// Set the counter upto the second heat (womb preparation (calving to first oestrous (Olds et al., 1953) + silent first heat)
									current_node->gest_at_inf= 0;
								} else {
									current_node->gest_at_inf= current_node->days_to_calv;
									current_node->calf_fate= _Fate_of_calf(current_node->inf_stat, current_node->gest_at_inf);
								}
							}
						}
					}
					previous_node= current_node;
					current_node= current_node->next_node;
				}
			}
		}
	}

}
