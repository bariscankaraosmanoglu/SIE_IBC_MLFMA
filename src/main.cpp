#include "main.h"

int main()
{
    read_simulationinput();
    generate_edgeconnection();
    generate_octtree();
    generate_nearfarboxes();
    generate_Znearcomponent();
    generate_Zfarcomponents();
    generate_translationcomponent();

    calculate_loadbalance();
    calculate_Znearmatrix();
    calculate_translationoperator();
    calculate_blockdiagonalpreconditioner();

    stage_solution();
    stage_aftersolution();
    clear_memory();



    return 0;
}
