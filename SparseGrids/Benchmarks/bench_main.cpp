#include "benchEvaluate.hpp"

void printHelp();

int main(int argc, char** argv){

    //cout << " Phruuuuphrrr " << endl; // this is the sound that the Tasmanian devil makes

    std::vector<std::string> args = stringArgs(argc, argv);

    if (args.empty() || hasHelp(args.back())){
        printHelp();
        return 0;
    }

    auto test = getTest(args.back());
    if (test == test_none){
        printHelp();
        return 1;
    }


    cout << "Testing eval: " << endl;

    return 0;
}

void printHelp(){
    cout << "\nusage: ./benchmark <function> <parameters>\n\n";
    cout << "functions: evaluate\n";
    cout << "\n see: ./benchmark <function> help\n\n";
}
