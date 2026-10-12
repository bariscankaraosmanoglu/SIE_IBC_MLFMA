// Copyright (c) 2026 Bariscan Karaosmanoglu, PhD
// SPDX-License-Identifier: MIT

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
    if (val_pc == 1) calculate_blockdiagonalpreconditioner();

    const int solution_status = stage_solution();
    stage_aftersolution();
    clear_memory();

    return solution_status;
}
