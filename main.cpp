#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <tuple>
#include <map>
#include <chrono>
#include <random>
#include <functional>
using namespace std;   

const size_t P = 600;
const double R = 0.167;
const double M = 0.33;
const bool VERBOSE = false;
const size_t POP_SIZE = P;
const size_t NUM_KILLED_PER_GEN = static_cast<int>(R * P);
const size_t NUM_MUTATED_PER_GEN = static_cast<int>(M * P);
const size_t SELECTOR_BUCKET_COUNT = NUM_MUTATED_PER_GEN + NUM_KILLED_PER_GEN;
const size_t MULTIPLIER = 1;
const size_t ITERATION_LIMIT = 100;
const size_t GENE_LOWER_LIMIT = 0;
const size_t GENE_UPPER_LIMIT = 75;

const float infinity = numeric_limits<float>::infinity();
unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine randGenerator(seed);
uniform_real_distribution<float> genUniform(0, 1);
uniform_real_distribution<float> genGene(GENE_LOWER_LIMIT, GENE_UPPER_LIMIT);


// Class:   Chromosome
// Members: int id, float fitness & 4 floats (w,w,y,z)
// Methods:
//          Chromosome() Default Constructor
//          Chromosome() Parameterized Constructor
//          print()
// Notes:   none;
struct Chromosome {
    int id = rand();
    float fitness;
    float w;
    float x;
    float y;
    float z;


    // Function:    Chromosome() Default Constructor
    // Parameters:  None
    // Returns:     Chromosome object
    // Purpose:     Initialize with silent Not A Number
    Chromosome() :
        w{nanf("")}, x{nanf("")}, y{nanf("")}, z{nanf("")}, fitness{nanf("")}{}


    // Function:    Chromosome() Parameterized Constructor
    // Parameters:  4 floats (x-z)
    // Returns:     Chromosome object
    // Purpose:     Initialize with values for genes
    Chromosome(float w_, float x_, float y_, float z_) :
        w{w_}, x{x_}, y{y_}, z{z_}, fitness{nanf("")} {}
};




// This needs context of a Chromosome first
typedef function<float (const Chromosome&)> FitnessFunc;


// Function:    mirrorIntoRange()
// Parameters:  1 float
// Returns:     float
// Purpose:     Keep the float within 0 <= x <= 75
float mirrorIntoRange(float gene){
    while(gene < GENE_LOWER_LIMIT || gene > GENE_UPPER_LIMIT){
        if(gene < GENE_LOWER_LIMIT)
            gene *= -1;
        if(gene > GENE_UPPER_LIMIT)
            gene = (GENE_UPPER_LIMIT * 2) - gene;
    }
    return gene;
}


// Function:    mutateGene()
// Parameters:  1 float
// Returns:     float
// Purpose:     Mutate the gene and make sure it's in the correct range
float mutateGene(float gene){
    float delta = MULTIPLIER * (0.5 - genUniform(randGenerator));
    return mirrorIntoRange(gene + delta);
}


// Function:    fitness1()
// Parameters:  1 Chromosome const reference
// Returns:     float
// Purpose:     Determine the fitness of the Chromosome
float fitness1(const Chromosome& c){
    const float omega = M_PI * 2;
    const float sigma = 150;
    return (
        pow(sin(omega * (c.w + c.y)), 2) *
        pow(sin(omega * (c.x + c.z)), 2) *
        exp(-(c.w + c.x + c.y + c.z) / sigma)
    );
}


// Function:    fitness2()
// Parameters:  1 Chromosome const reference
// Returns:     float
// Purpose:     Determine the fitness of the Chromosome (alternate)
float fitness2(const Chromosome& c) {
    return static_cast<float>(1) - fitness1(c);
}


// Function:    findHighestRankingChromosome()
// Parameters:  1 vector of Chromosomes const references and a comparison
//              function
// Returns:     Chromosome
// Purpose:     Determine the Chromosome that ranks higher than all others,
//              based on the comparison function. This functions is only a
//              helper function to help create other functions.
Chromosome findHighestRankingChromosome(
    const vector<Chromosome>& population,
    function<bool(float,float)> compare
){
    Chromosome best = population[0];
    for(const Chromosome& c : population)
        if (compare(c.fitness, best.fitness))
            best = c;
    return best;
}


// Function:    findBestChromosome()
// Parameters:  1 vector of Chromosomes const references
// Returns:     Chromosome
// Purpose:     Determine the best Chromosome
Chromosome findBestChromosome(const vector<Chromosome>& population){
    return findHighestRankingChromosome(population, greater<float>());
}


// Function:    findWorstChromosome()
// Parameters:  1 vector of Chromosomes const references
// Returns:     Chromosome
// Purpose:     Determine the worst Chromosome
Chromosome findWorstChromosome(const vector<Chromosome>& population){
    return findHighestRankingChromosome(population, less<float>());
}


// Function:    sumTotalFitness()
// Parameters:  1 vector of Chromosomes const reference
// Returns:     float
// Purpose:     Calculate the total fitness for a generation
float sumTotalFitness(const vector<Chromosome>& chromosomes){
    float total = 0;
    for(const Chromosome& c : chromosomes)
        total += c.fitness;
    return total;
}




// Class:   Selector
// Members: 1 vector of maps of float and Chromosome
// Methods:
//          constructor()
//          selectOne()
// Notes:   none
class Selector {
    vector<map<float, const Chromosome&>> buckets;

public:
    // Function:    Selector() Constructor
    // Parameters:  1 vector of Chromosomes const reference & 1 const fitness
    //              function
    // Returns:     Selector object
    // Purpose:     Select the "fittest" of the generation
    Selector(const vector<Chromosome>& population){
        map<float, const Chromosome&> bucket;
        float currentFitnessPercent = 0;
        float totalFitness = sumTotalFitness(population);
        for(const Chromosome& c : population){
            currentFitnessPercent += c.fitness / totalFitness;
            while(
                currentFitnessPercent
                >=
                (
                    static_cast<float>(1 + buckets.size())
                    /
                    static_cast<float>(SELECTOR_BUCKET_COUNT)
                )
            ){
                bucket.insert(
                    pair<float, const Chromosome&>(currentFitnessPercent, c)
                );
                buckets.push_back(bucket);
                bucket = map<float, const Chromosome&>();
            }
            bucket.insert(
                pair<float, const Chromosome&>(currentFitnessPercent, c)
            );
        }
        while(buckets.size() < SELECTOR_BUCKET_COUNT){
            // This means the bucket isn't full due to floating point errors
            bucket.insert(pair<float, const Chromosome&>(1, population[0]));
            // Adds chromosome to fill the gap created by floating point error.
            buckets.push_back(bucket);
            bucket = map<float, const Chromosome&>();
        }
    }


    // Function:    selectOne()
    // Parameters:  None
    // Returns:     Chromosome
    // Purpose:     Select one Chromosome from the generation
    Chromosome selectOne() const{
        float n = genUniform(randGenerator);
        const map<float, const Chromosome&>& bucket = buckets[
            size_t(n * SELECTOR_BUCKET_COUNT)
        ];
        // This is decomposition to get key-value pair
        // (C++17 feature)
        for(auto const& [limit, c] : bucket){
            if (n <= limit)
                return c;
        }
        // There's no catch; this should cause the program to fail
        throw runtime_error("Unreachable");
    }
};




// Function:    selectSurvivors()
// Parameters:  vector of Chromosome const references
// Returns:     vector of Chromosomes
// Purpose:     Select the "survivors" or "best" of a generation
vector<Chromosome> selectSurvivors(const vector<Chromosome>& population){
    Selector selector{population};
    vector<Chromosome> result;
    for(size_t i = 0; i < (POP_SIZE - NUM_KILLED_PER_GEN); ++i)
        result.push_back(selector.selectOne());
    return result;
}


// Function:    crossover()
// Parameters:  2 Chromosome const references & fitness function
// Returns:     Chromosome
// Purpose:     Crossover the genes of the Chromosome parents to make child
Chromosome crossover(
    const Chromosome& mate1,
    const Chromosome& mate2,
    const FitnessFunc fitness
){
    float w = genUniform(randGenerator) > 0.5 ? mate1.w : mate2.w;
    float x = genUniform(randGenerator) > 0.5 ? mate1.x : mate2.x;
    float y = genUniform(randGenerator) > 0.5 ? mate1.y : mate2.y;
    float z = genUniform(randGenerator) > 0.5 ? mate1.z : mate2.z;
    Chromosome child = Chromosome(w, x, y, z);
    child.fitness = fitness(child);
    return child;
}


// Function:    reproduce()
// Parameters:  vector of Chromosome const references & fitness function
// Returns:     vector of Chromosomes
// Purpose:     Reproduce and replace a portion of a generation
vector<Chromosome> reproduce(
    const vector<Chromosome>& population,
    const FitnessFunc fitness
){
    Selector selector{population};
    vector<Chromosome> result{population};
    while(result.size() < POP_SIZE){
        Chromosome mate1 = selector.selectOne();
        Chromosome mate2 = selector.selectOne();
        for(size_t i = 0; i < 2; ++i){
            Chromosome child = crossover(mate1, mate2, fitness);
            result.push_back(child);
        }
    }
    return result;
}





// Function:    mutatePopulation()
// Parameters:  vector of Chromosome const references & fitness function
// Returns:     vector of Chromosomes
// Purpose:     Mutate a specified number of Chromosomes per generation and
//              return all of those Chromosomes
vector<Chromosome> mutatePopulation(
    const vector<Chromosome>& population,
    const FitnessFunc fitness
){
    vector<Chromosome> result{population};
    vector<Chromosome> toMutate{};
    for(size_t i = 0; i < NUM_MUTATED_PER_GEN; ++i){
        size_t randIndex = static_cast<size_t>(
            result.size() * genUniform(randGenerator)
        );
        Chromosome element = *next(result.begin(), randIndex);
        result.erase(next(result.begin(), randIndex));
        toMutate.push_back(element);
    }
    for(Chromosome& c : toMutate){
        Chromosome child = Chromosome(
            mutateGene(c.w),
            mutateGene(c.x),
            mutateGene(c.y),
            mutateGene(c.z)
        );
        child.fitness = fitness(child);
        result.push_back(child);
    }
    return result;
}


// Function:    nextGeneration()
// Parameters:  vector of Chromosome const references & fitness function
// Returns:     vector of Chromosomes
// Purpose:     Create a new generation of Chromosomes by first selecting
//              survivors and then mutating a percentage of the population
vector<Chromosome> nextGeneration(
    const vector<Chromosome>& population_,
    const FitnessFunc fitness
){
    vector<Chromosome> population = selectSurvivors(population_);
    population = reproduce(population, fitness);
    return mutatePopulation(population, fitness);
}


// Function:    generateRandomPopulation()
// Parameters:  None
// Returns:     vector of Chromosomes
// Purpose:     Create a seed generation of Chromosomes
vector<Chromosome> generateRandomPopulation(const FitnessFunc fitness){
    vector<Chromosome> population{};
    for(size_t i = 0; i < POP_SIZE; ++i){
        Chromosome child = Chromosome(
            genGene(randGenerator),
            genGene(randGenerator),
            genGene(randGenerator),
            genGene(randGenerator)
        );
        child.fitness = fitness(child);
        population.push_back(child);
    }
    return population;
}


// Function:    run()
// Parameters:  vector of Chromosome reference & fitness function
// Returns:     The population containing the best chromosome
// Purpose:     Determine the top Chromosome of a generation
vector<Chromosome> run(
    vector<Chromosome>& population,
    const FitnessFunc fitness
){
    vector<Chromosome> topPop;
    float topFitness = -numeric_limits<float>::infinity();
    for(size_t i = 0; i < ITERATION_LIMIT; ++i){
        if (i>0)
            population = nextGeneration(population, fitness);
        Chromosome maybeBest = findBestChromosome(population);
        if(maybeBest.fitness > topFitness){
            topPop = population;
            topFitness = maybeBest.fitness;
        }
    }
    return topPop;
}


// Function:    print()
// Parameters:  vector of Chromosome references & size_t iteration
// Returns:     void
// Purpose:     print out the genetic info
void print(const vector<Chromosome>& pool, size_t iteration){
    cout << "\nBest generation from run #" << iteration << endl;
    cout << setw(10) << "Chromosome" << "  ";
    cout << setw(10) << "Gene W" << "   ";
    cout << setw(10) << "Gene X" << "   ";
    cout << setw(10) << "Gene Y" << "   ";
    cout << setw(10) << "Gene Z" << "  ";
    cout << setw(11) << "Fitness" << endl;
    for(size_t j = 0; j < pool.size(); ++j){
        cout << setw(10) << j + 1 << ") ";
        cout << setw(10) << pool[j].w << " | ";
        cout << setw(10) << pool[j].x << " | ";
        cout << setw(10) << pool[j].y << " | ";
        cout << setw(10) << pool[j].z << " | ";
        cout << setw(11) << pool[j].fitness << endl;
    }
}


// Function:    printPopulationStatistics()
// Parameters:  vector of Chromosome references
// Returns:     void
// Purpose:     Prints statistics about a population
void printPopulationStatistics(const vector<Chromosome>& pool){
    float average = static_cast<float>(sumTotalFitness(pool) / pool.size());
    float stddev = 0;
    for(const Chromosome c : pool)
        stddev += pow(average - c.fitness, 2);
    stddev = pow((stddev / pool.size()), .5);
    cout << "best: " << findBestChromosome(pool).fitness;
    cout << " worst: " << findWorstChromosome(pool).fitness;
    cout << " average: " << average;
    cout << " stdev: " << stddev << endl;
}


// Function:    test()
// Parameters:  vector of Chromosome references & fitness function & optional
//              count (of generations to go through)
// Returns:     void
// Purpose:     test out the genetic algorithm with a fitness function and
//              maintain the results
void test(
    vector<Chromosome>& results_,
    const FitnessFunc fitness,
    const string& name,
    size_t count,
    bool debug=false
){
    results_.clear();
    cout << "\nTesting with " << name << endl;
    for(size_t i = 0; i < count; ++i){
        vector<Chromosome> population = generateRandomPopulation(fitness);
        vector<Chromosome> topPop = run(population, fitness);
        Chromosome top = findBestChromosome(topPop);
        if(debug)
            print(topPop, (i + 1));
        cout << "Statistics from the best generation in " << i + 1 << ":";
        cout <<  endl;
        printPopulationStatistics(topPop);
        cout << "Best chromosomes: w+y=" << top.w + top.y
             << " x+z=" << top.x + top.z << endl << endl;
        results_.push_back(top);
    }
    cout << "---------- Statistics across all runs with this fitness function";
    cout << " ----------" << endl;
    printPopulationStatistics(results_);
    cout << endl;
}





//---------------------------- Main Method ----------------------------------//
int main(){
    size_t count = 10;
    vector<Chromosome> results;
    test(results, fitness1, "provided fitness function", count, VERBOSE);
    test(results, fitness2, "modified fitness function", count, VERBOSE);
    return 0;
}
