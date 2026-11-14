/*
 * babymodel custom: monomer -> oligomer kinetics, plaque nucleation/growth, coloring
 */

#include "./custom.h"
#include <algorithm>
#include <cmath>

static int substrate_index( const char* name )
{
	int idx = microenvironment.find_density_index( name );
	if( idx < 0 )
	{ std::cout << "Warning: substrate " << name << " not found." << std::endl; }
	return idx;
}

static void update_monomer_at_neurons( double dt )
{
	static int monomer_idx = -1;
	static int oligomer_idx = -1;
	static Cell_Definition* pNeuron = NULL;
	if( monomer_idx < 0 )
	{
		monomer_idx = substrate_index("amyloid_monomer");
		oligomer_idx = substrate_index("amyloid_oligomer");
		pNeuron = find_cell_definition("neuron");
	}
	if( !pNeuron || monomer_idx < 0 || oligomer_idx < 0 )
	{ return; }

	double threshold = 0.04;
	double rate = 0.45;
	if( parameters.doubles.find_index("monomer_oligo_threshold") != -1 )
	{ threshold = parameters.doubles("monomer_oligo_threshold"); }
	if( parameters.doubles.find_index("monomer_oligo_rate") != -1 )
	{ rate = 2.0 * parameters.doubles("monomer_oligo_rate"); }

	for( auto* pC : *all_cells )
	{
		if( pC->type != pNeuron->type || pC->phenotype.death.dead )
		{ continue; }
		int v = microenvironment.nearest_voxel_index( pC->position );
		if( v < 0 )
		{ continue; }
		double& m = microenvironment.density_vector(v)[monomer_idx];
		if( m <= threshold )
		{ continue; }
		double excess = m - threshold;
		double transfer = rate * dt * excess;
		transfer = std::min( transfer, m * 0.5 );
		microenvironment.density_vector(v)[oligomer_idx] += transfer;
		m -= transfer;
	}
}

static void update_monomer_to_oligomer( double dt )
{
	static int monomer_idx = -1;
	static int oligomer_idx = -1;
	if( monomer_idx < 0 )
	{
		monomer_idx = substrate_index("amyloid_monomer");
		oligomer_idx = substrate_index("amyloid_oligomer");
	}
	if( monomer_idx < 0 || oligomer_idx < 0 )
	{ return; }

	double threshold = 0.04;
	double rate = 0.18;
	if( parameters.doubles.find_index("monomer_oligo_threshold") != -1 )
	{ threshold = parameters.doubles("monomer_oligo_threshold"); }
	if( parameters.doubles.find_index("monomer_oligo_rate") != -1 )
	{ rate = parameters.doubles("monomer_oligo_rate"); }

	#pragma omp parallel for
	for( int n = 0; n < microenvironment.number_of_voxels(); n++ )
	{
		double m = microenvironment.density_vector(n)[monomer_idx];
		if( m <= threshold )
		{ continue; }
		double excess = m - threshold;
		double transfer = rate * dt * excess;
		transfer = std::min( transfer, m * 0.35 );
		transfer = std::min( transfer, excess * 0.5 );
		microenvironment.density_vector(n)[oligomer_idx] += transfer;
		microenvironment.density_vector(n)[monomer_idx] -= transfer;
	}
}

static void subtract_at_position( std::vector<double>& position, int substrate_idx, double amount )
{
	if( substrate_idx < 0 || amount <= 0.0 )
	{ return; }
	int v = microenvironment.nearest_voxel_index( position );
	if( v < 0 )
	{ return; }
	double& dens = microenvironment.density_vector(v)[substrate_idx];
	dens = std::max( 0.0, dens - amount );
}

static std::vector< std::vector<double> > pending_plaques;

static void try_plaque_nucleation( double dt )
{
	static int oligomer_idx = -1;
	static int monomer_idx = -1;
	static Cell_Definition* pNeuron = NULL;
	static Cell_Definition* pPlaque = NULL;
	if( oligomer_idx < 0 )
	{
		oligomer_idx = substrate_index("amyloid_oligomer");
		monomer_idx = substrate_index("amyloid_monomer");
		pNeuron = find_cell_definition("neuron");
		pPlaque = find_cell_definition("amyloid_plaque");
	}
	if( !pNeuron || !pPlaque || oligomer_idx < 0 )
	{ return; }

	double o_thresh = 0.038;
	double hazard = 0.055;
	if( parameters.doubles.find_index("plaque_oligo_threshold") != -1 )
	{ o_thresh = parameters.doubles("plaque_oligo_threshold"); }
	if( parameters.doubles.find_index("plaque_nucleation_prob") != -1 )
	{ hazard = parameters.doubles("plaque_nucleation_prob"); }

	for( auto* pC : *all_cells )
	{
		if( pC->type != pNeuron->type || pC->phenotype.death.dead )
		{ continue; }

		int v = microenvironment.nearest_voxel_index( pC->position );
		if( v < 0 )
		{ continue; }
		double o = microenvironment.density_vector(v)[oligomer_idx];
		double m = microenvironment.density_vector(v)[monomer_idx];

		// One plaque per neuron neighborhood. Without this, every neuron
		// above threshold seeds a new plaque every few minutes.
		bool plaque_nearby = false;
		for( unsigned int k = 0; k < all_cells->size(); k++ )
		{
			Cell* pP = (*all_cells)[k];
			if( pP->type != pPlaque->type || pP->phenotype.death.dead )
			{ continue; }
			double dx = pP->position[0] - pC->position[0];
			double dy = pP->position[1] - pC->position[1];
			if( dx*dx + dy*dy < 40.0 * 40.0 )
			{
				plaque_nearby = true;
				break;
			}
		}
		for( unsigned int k = 0; k < pending_plaques.size() && !plaque_nearby; k++ )
		{
			double dx = pending_plaques[k][0] - pC->position[0];
			double dy = pending_plaques[k][1] - pC->position[1];
			if( dx*dx + dy*dy < 40.0 * 40.0 )
			{ plaque_nearby = true; }
		}
		if( plaque_nearby )
		{ continue; }

		bool o_ok = ( o >= o_thresh );
		bool m_ok = ( m >= 0.075 && o >= 0.02 );
		if( !o_ok && !m_ok )
		{ continue; }
		// hazard [1/min] * dt [min]; scale by local oligomer above threshold
		double o_excess = std::max( 0.0, o - o_thresh );
		double p_event = hazard * dt * ( 1.0 + 4.0 * o_excess );
		if( m_ok && !o_ok )
		{ p_event *= 0.65; }
		if( UniformRandom() > p_event )
		{ continue; }

		std::vector<double> pos = pC->position;
		pos[0] += 10.0 * ( UniformRandom() - 0.5 );
		pos[1] += 10.0 * ( UniformRandom() - 0.5 );

		// Queue the agent. create_cell() inside the parallel cell loop
		// reallocates all_cells and crashes once plaques start appearing.
		pending_plaques.push_back( pos );

		double mass = 0.05 * o;
		subtract_at_position( pos, oligomer_idx, mass );
		subtract_at_position( pos, monomer_idx, 0.5 * mass );
	}
}

void create_cell_types( void )
{
	if (parameters.ints.find_index("random_seed") != -1)
	{ SeedRandom(parameters.ints("random_seed")); }

	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	cell_defaults.functions.update_migration_bias = NULL;
	cell_defaults.functions.update_phenotype = NULL;
	cell_defaults.functions.custom_cell_rule = NULL;
	cell_defaults.functions.contact_function = NULL;
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	initialize_cell_definitions_from_pugixml();
	build_cell_definitions_maps();
	setup_signal_behavior_dictionaries();
	setup_cell_rules();

	cell_defaults.functions.update_phenotype = phenotype_function;
	cell_defaults.functions.custom_cell_rule = custom_function;
	cell_defaults.functions.contact_function = contact_function;

	// XML-loaded cell definitions keep their own function pointers; copy custom hooks to every type.
	for( int n = 0; n < cell_definitions_by_index.size(); n++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[n];
		if( pCD == NULL )
		{ continue; }
		pCD->functions.update_phenotype = phenotype_function;
		pCD->functions.custom_cell_rule = custom_function;
		pCD->functions.contact_function = contact_function;
	}

	display_cell_definitions( std::cout );
	return;
}

void setup_microenvironment( void )
{
	initialize_microenvironment();
	return;
}

void setup_tissue( void )
{
	load_cells_from_pugixml();
	set_parameters_from_distributions();

	std::mt19937 rng(42);
	std::uniform_real_distribution<double> plaque_vol_dist(2000.0, 12000.0);

	Cell_Definition* pPlaqueDef = find_cell_definition("amyloid_plaque");
	Cell_Definition* pSenescentDef = find_cell_definition("senescent_microglia");
	Cell_Definition* pNeuronDef = find_cell_definition("neuron");

	for( auto* pC : *all_cells )
	{
		if( pPlaqueDef != nullptr && pC->type == pPlaqueDef->type )
		{
			pC->set_total_volume( plaque_vol_dist(rng) );
			pC->is_movable = false;
		}
		if( pSenescentDef != nullptr && pC->type == pSenescentDef->type )
		{ pC->is_movable = false; }
		if( pNeuronDef != nullptr && pC->type == pNeuronDef->type )
		{ pC->is_movable = false; }
	}
	return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector<std::string> output = { "black", "black", "black", "black" };

	if( pCell->phenotype.death.dead )
	{
		output[0] = "rgb(200,200,200)";
		output[2] = "rgb(100,100,100)";
		return output;
	}

	static Cell_Definition* pInact = NULL;
	static Cell_Definition* pAct = NULL;
	static Cell_Definition* pPlaque = NULL;
	static Cell_Definition* pSen = NULL;
	static Cell_Definition* pNeuron = NULL;
	if( !pInact )
	{
		pInact = find_cell_definition("inactivated_microglia");
		pAct = find_cell_definition("activated_microglia");
		pPlaque = find_cell_definition("amyloid_plaque");
		pSen = find_cell_definition("senescent_microglia");
		pNeuron = find_cell_definition("neuron");
	}

	if( pInact && pCell->type == pInact->type )
	{
		output[0] = "rgb(100,149,237)";
		output[2] = "rgb(30,60,120)";
	}
	else if( pAct && pCell->type == pAct->type )
	{
		output[0] = "rgb(220,60,60)";
		output[2] = "rgb(120,20,20)";
	}
	else if( pPlaque && pCell->type == pPlaque->type )
	{
		double frac = std::min(1.0, (pCell->phenotype.volume.total - 2000.0) / 10000.0);
		int r = (int)(240 - frac * 80);
		int g = (int)(190 - frac * 100);
		int b = (int)(20 + frac * 5);
		output[0] = "rgb(" + std::to_string(r) + "," + std::to_string(g) + "," + std::to_string(b) + ")";
		output[2] = "rgb(" + std::to_string(r/2) + "," + std::to_string(g/2) + "," + std::to_string(b/2) + ")";
	}
	else if( pSen && pCell->type == pSen->type )
	{
		output[0] = "rgb(140,80,160)";
		output[2] = "rgb(70,30,90)";
	}
	else if( pNeuron && pCell->type == pNeuron->type )
	{
		output[0] = "rgb(60,180,80)";
		output[2] = "rgb(20,90,40)";
	}

	return output;
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	static Cell_Definition* pPlaque = NULL;
	static Cell_Definition* pSen = NULL;
	static Cell_Definition* pInact = NULL;
	static Cell_Definition* pAct = NULL;
	static int oligomer_idx = -1;
	if( !pPlaque )
	{
		pPlaque = find_cell_definition("amyloid_plaque");
		pSen = find_cell_definition("senescent_microglia");
		pInact = find_cell_definition("inactivated_microglia");
		pAct = find_cell_definition("activated_microglia");
		oligomer_idx = substrate_index("amyloid_oligomer");
	}

	if( pSen && pCell->type == pSen->type )
	{ pCell->is_movable = false; }
	if( pPlaque && pCell->type == pPlaque->type )
	{ pCell->is_movable = false; }

	if( pPlaque && pCell->type == pPlaque->type )
	{
		double frac = (pCell->phenotype.volume.total - 2000.0) / 10000.0;
		frac = std::max(0.0, std::min(1.0, frac));

		// Plaques do not secrete Aβ; insoluble deposit only.
		phenotype.secretion.sync_to_microenvironment( &microenvironment );
		int monomer_idx = substrate_index("amyloid_monomer");
		if( monomer_idx >= 0 )
		{ phenotype.secretion.secretion_rates[monomer_idx] = 0.0; }
		if( oligomer_idx >= 0 )
		{ phenotype.secretion.secretion_rates[oligomer_idx] = 0.0; }

		phenotype.mechanics.cell_cell_repulsion_strength = 2.0 + frac * 6.0;
		phenotype.mechanics.cell_cell_adhesion_strength = 0.1 + frac * 0.5;

		if( pInact )
		{ phenotype.mechanics.cell_adhesion_affinities[pInact->type] = 0.2 + frac * 0.4; }
		if( pAct )
		{ phenotype.mechanics.cell_adhesion_affinities[pAct->type] = 1.0 + frac * 2.5; }
		if( pPlaque )
		{ phenotype.mechanics.cell_adhesion_affinities[pPlaque->type] = 0.3 + frac * 1.8; }
		if( pSen )
		{ phenotype.mechanics.cell_adhesion_affinities[pSen->type] = 0.0; }

		// Grow plaque by consuming local oligomer (reaccumulation after partial clearance).
		double growth = 0.0004;
		if( parameters.doubles.find_index("plaque_growth_rate") != -1 )
		{ growth = parameters.doubles("plaque_growth_rate"); }

		int v = microenvironment.nearest_voxel_index( pCell->position );
		if( v >= 0 && oligomer_idx >= 0 )
		{
			double o = microenvironment.density_vector(v)[oligomer_idx];
			double uptake = growth * dt * o;
			uptake = std::min( uptake, o * 0.1 );
			microenvironment.density_vector(v)[oligomer_idx] -= uptake;
			double new_vol = pCell->phenotype.volume.total + uptake * 5000.0;
			new_vol = std::min( new_vol, 12000.0 );
			pCell->set_total_volume( new_vol );
		}
	}
	return;
}

void spawn_queued_plaques( void )
{
	if( pending_plaques.empty() )
	{ return; }
	Cell_Definition* pPlaque = find_cell_definition( "amyloid_plaque" );
	if( pPlaque == NULL )
	{
		pending_plaques.clear();
		return;
	}
	for( unsigned int n = 0; n < pending_plaques.size(); n++ )
	{
		Cell* pNew = create_cell( *pPlaque );
		pNew->assign_position( pending_plaques[n] );
		pNew->set_total_volume( 2500.0 );
		pNew->is_movable = false;
	}
	pending_plaques.clear();
	return;
}

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
	// Global kinetics once per mechanics step (OpenMP-safe).
	bool run_kinetics = false;
	#pragma omp critical( babymodel_global_kinetics )
	{
		static double last_mechanics_time = -1.0;
		if( PhysiCell_globals.current_time != last_mechanics_time )
		{
			last_mechanics_time = PhysiCell_globals.current_time;
			run_kinetics = true;
		}
	}
	if( run_kinetics )
	{
		update_monomer_to_oligomer( dt );
		update_monomer_at_neurons( dt );
		try_plaque_nucleation( dt );

		// Progress log every 120 simulated minutes (visible in Studio Run console).
		static double last_log_time = -1.0;
		double t = PhysiCell_globals.current_time;
		if( t - last_log_time >= 1440.0 || last_log_time < 0.0 )
		{
			last_log_time = t;
			int n_plaque = 0;
			int n_sen = 0;
			static Cell_Definition* pPlaque = NULL;
			static Cell_Definition* pSen = NULL;
			if( !pPlaque )
			{
				pPlaque = find_cell_definition("amyloid_plaque");
				pSen = find_cell_definition("senescent_microglia");
			}
			for( auto* pC : *all_cells )
			{
				if( pPlaque && pC->type == pPlaque->type && !pC->phenotype.death.dead )
				{ n_plaque++; }
				if( pSen && pC->type == pSen->type && !pC->phenotype.death.dead )
				{ n_sen++; }
			}
			double max_o = 0.0;
			double max_m = 0.0;
			int o_idx = substrate_index("amyloid_oligomer");
			int m_idx = substrate_index("amyloid_monomer");
			for( int k = 0; k < microenvironment.number_of_voxels(); k++ )
			{
				if( m_idx >= 0 )
				{ max_m = std::max( max_m, microenvironment.density_vector(k)[m_idx] ); }
				if( o_idx >= 0 )
				{ max_o = std::max( max_o, microenvironment.density_vector(k)[o_idx] ); }
			}
			std::cout << "[babymodel t=" << t << " min] amyloid_plaque=" << n_plaque
			          << " senescent_microglia=" << n_sen
			          << " total_agents=" << all_cells->size()
			          << " max_monomer=" << max_m << " max_oligomer=" << max_o << std::endl;
		}
	}
	return;
}

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{
	static Cell_Definition* pAct = NULL;
	static Cell_Definition* pPlaque = NULL;
	if( !pAct )
	{
		pAct = find_cell_definition("activated_microglia");
		pPlaque = find_cell_definition("amyloid_plaque");
	}
	if( !pAct || !pPlaque )
	{ return; }

	// Activated microglia erode plaque volume on sustained contact (complement to phagocytosis).
	if( pMe->type == pAct->type && pOther->type == pPlaque->type && !pOther->phenotype.death.dead )
	{
		double shrink = 0.45 * dt;
		double new_vol = std::max( 1800.0, pOther->phenotype.volume.total - shrink );
		pOther->set_total_volume( new_vol );
	}
	else if( pOther->type == pAct->type && pMe->type == pPlaque->type && !pMe->phenotype.death.dead )
	{
		double shrink = 0.45 * dt;
		double new_vol = std::max( 1800.0, pMe->phenotype.volume.total - shrink );
		pMe->set_total_volume( new_vol );
	}
	return;
}
