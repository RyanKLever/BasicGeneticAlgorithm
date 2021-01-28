#include <cmath>
#include <iostream>
#include <vector>
#include <list>
#include <tuple>
#include <map>
#include <chrono>
#include <random>
#include <functional>
using namespace std;

const int POP_SIZE = 600;
const int NUM_KILLED_PER_GEN = 100;
const int NUM_MUTATED_PER_GEN = 200;
const float STARTING_MUTATION_SCALE = 50;
const float MUTATION_SCALE_FACTOR = 0.985;
const float END_THRESHOLD = 0.03;
const float WHEN_NORMALIZING_IGNORE_LOWEST = NUM_MUTATED_PER_GEN + NUM_KILLED_PER_GEN;
const int SELECTOR_BUCKET_COUNT = POP_SIZE/2;

const float infinity = numeric_limits<float>::infinity();

struct Chromosome {
    int id = rand();
    float w;
    float x;
    float y;
    float z;
    Chromosome() : w{nanf("")}, x{nanf("")}, y{nanf("")}, z{nanf("")} {}
    Chromosome(float w_, float x_, float y_, float z_) :
    w{w_}, x{x_}, y{y_}, z{z_} {}
};

typedef function<float (const Chromosome&)> FitnessFunc;

unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine randGenerator(seed);
uniform_real_distribution<float> genUniform(0, 1);
uniform_real_distribution<float> genGene(0, 75);

float mirrorIntoRange(float gene){
  while(gene < 0 || gene > 75) {
    if (gene < 0)
      gene *= -1;
    if(gene > 75)
      gene = 150 - gene;
  }
  return gene;
}

float mutateGene(float gene, float multiplier){
    float delta = multiplier * (0.5 - genUniform(randGenerator));
    return mirrorIntoRange(gene + delta);
}

float fitness1(const Chromosome& c){
    const float omega = M_PI * 2;
    const float sigma = 150;
    return (
        pow(sin(omega * (c.w + c.y)), 2) *
        pow(sin(omega * (c.x + c.z)), 2) *
        exp(-(c.w + c.x + c.y + c.z) / sigma)
    );
}

float fitness2(const Chromosome& c){
    return 1 - fitness1(c);
}

Chromosome findBestChromosome(
    const vector<Chromosome>& population,
    FitnessFunc fitness
){
    Chromosome best = population[0];
    for (const Chromosome& c : population)
        if (fitness(best) < fitness(c))
            best = c;
    return best;
}

float sumTotalFitness(
    const vector<Chromosome>& chromosomes,
    FitnessFunc fitness
){
    float total = 0.0;
    for (const Chromosome& c : chromosomes)
        total += fitness(c);
    return total;
}

class Selector {
    vector<map<float, const Chromosome&>> buckets;

public:
    Selector(
        const vector<Chromosome>& population,
        const FitnessFunc normalizedFitness
    ){
        map<float, const Chromosome&> bucket;
        float currentFitnessPercent = 0.0;
        float totalFitness = sumTotalFitness(population, normalizedFitness);
        for (const Chromosome& c : population) {
            if (normalizedFitness(c) == 0.0)
                continue;
            currentFitnessPercent += normalizedFitness(c) / totalFitness;
            while (
                currentFitnessPercent
                >=
                (1+buckets.size()) / float(SELECTOR_BUCKET_COUNT)
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
        while (buckets.size() < SELECTOR_BUCKET_COUNT){
            // The last bucket wasn't incuded due to floating point errors
            bucket.insert(pair<float, const Chromosome&>(1, population[0]));
            // Adds a chromosome to fill loss created by floating point error.
            buckets.push_back(bucket);
            bucket = map<float, const Chromosome&>();
        }
    }

    Chromosome selectOne() const {
        float n = genUniform(randGenerator);
        const map<float, const Chromosome&>& bucket = buckets[
            size_t(n * SELECTOR_BUCKET_COUNT)
        ];

        for (auto const& [limit, c] : bucket)
            if (n <= limit)
                return c;

        throw runtime_error("unreachable");
    }
};

FitnessFunc normalizeFitnessFunc(
    const FitnessFunc fitness,
    const vector<Chromosome>& population
){
    vector<float> fitnesses;
    for (const Chromosome& c: population)
        fitnesses.push_back(fitness(c));
    sort(fitnesses.begin(), fitnesses.end());

    float highest = fitnesses[fitnesses.size()-1];
    float lowest = fitnesses[WHEN_NORMALIZING_IGNORE_LOWEST];
    // We exclude the most lowest points
    if (highest == lowest)
        return [](const Chromosome& c) -> float { return 1.0; };

    return [highest, lowest, &fitness](const Chromosome& c) -> float {
        float score = (fitness(c)-lowest) / (highest-lowest);
        return score < 0 ? 0 : score;
    };
}

vector<Chromosome> selectSurvivors(
    const vector<Chromosome>& population,
    const FitnessFunc normalizedFitness
){
    Selector selector{population, normalizedFitness};

    vector<Chromosome> result;
    for (int i=0; i<POP_SIZE-NUM_KILLED_PER_GEN; ++i)
        result.push_back( selector.selectOne() );
    return result;
}

vector<Chromosome> repopulateByCrossover(
    const vector<Chromosome>& population,
    const FitnessFunc normalizedFitness
){
    Selector selector{population, normalizedFitness};

    vector<Chromosome> result{population};
    while (result.size() < POP_SIZE) {
        Chromosome mate1 = selector.selectOne();
        Chromosome mate2 = selector.selectOne();
        for (int i = 0; i < 2; ++i) {
            float w = genUniform(randGenerator) > 0.5 ? mate1.w : mate2.w;
            float x = genUniform(randGenerator) > 0.5 ? mate1.x : mate2.x;
            float y = genUniform(randGenerator) > 0.5 ? mate1.y : mate2.y;
            float z = genUniform(randGenerator) > 0.5 ? mate1.z : mate2.z;
            result.push_back(Chromosome(w, x, y, z) );
        }
    }

    return result;
}

vector<Chromosome> mutatePopulation(
    const vector<Chromosome>& population,
    float mutatationScale
){
    vector<Chromosome> result{population};
    vector<Chromosome> toMutate{};

    for(int i = 0; i < NUM_MUTATED_PER_GEN; ++i){
        size_t randIndex = int(result.size() * genUniform(randGenerator));
        Chromosome element = *next(result.begin(), randIndex);
        result.erase(next(result.begin(), randIndex));
        toMutate.push_back(element);
    }

    for(Chromosome& c : toMutate){
        result.push_back(
            Chromosome(
                mutateGene(c.w, mutatationScale),
                mutateGene(c.x, mutatationScale),
                mutateGene(c.y, mutatationScale),
                mutateGene(c.z, mutatationScale)
            )
        );
    }

    return result;
}

vector<Chromosome> nextGeneration(
    const vector<Chromosome>& population_,
    const FitnessFunc fitness,
    float mutationScale
){
    FitnessFunc normalizedFitness = normalizeFitnessFunc(fitness, population_);

    vector<Chromosome> population = selectSurvivors(
        population_,
        normalizedFitness
    );
    population = repopulateByCrossover(population, normalizedFitness);
    return mutatePopulation(population, mutationScale);
}

vector<Chromosome> generateRandomPopulation(){
    vector<Chromosome> population{};
    for(int i = 0; i < POP_SIZE; ++i)
        population.push_back(
            Chromosome(
                genGene(randGenerator), genGene(randGenerator),
                genGene(randGenerator), genGene(randGenerator)
            )
        );
    return population;
}

Chromosome run(
    const vector<Chromosome>& population_,
    const FitnessFunc fitness
){
    vector<Chromosome> population{population_};
    float mutationScale = STARTING_MUTATION_SCALE;
    Chromosome top;
    float topFitness = -numeric_limits<float>::infinity();
    while(true){
        Chromosome maybeBest = findBestChromosome(population, fitness);
        if(fitness(maybeBest) > topFitness){
            top = maybeBest;
            topFitness = fitness(maybeBest);
        }

        mutationScale *= MUTATION_SCALE_FACTOR;
        if (mutationScale <= END_THRESHOLD)
            return top;

        population = nextGeneration(population, fitness, mutationScale);
    }
}

int main() {
    int count = 10;
    vector<Chromosome> results;
    for (int i=0; i<count; ++i) {
        vector<Chromosome> population = generateRandomPopulation();
        Chromosome top = run(population, fitness1);
        cout << fitness1(top) << " w+y=" << top.w + top.y << " x+z=";
        cout << top.x + top.z << endl;
        results.push_back(top);
    }

    float average = sumTotalFitness(results, fitness1) / count;
    float stddev = 0;
    for (Chromosome c : results)
        stddev += pow(average - fitness1(c), 2);
    stddev = pow(stddev/count, .5);
    cout << "average: " << average << " stdev: " << stddev << endl;
    return 0;
}
