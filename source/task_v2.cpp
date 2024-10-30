// Input/output
#include <stdio.h>
#include <iostream>
#include <fstream>
// sin, cos...
#include <cmath>
// Strings
#include <string>
#include <cstring>
// Time
#include <chrono>
// MPI
#include <mpi.h>

using namespace std;


// Check if point is in D
bool is_point_in_D(double x, double y)
{
    // Triangle is given by 3 lines, thus the given point should be:
    // 1) at or above the line from (3,0) to (3,0)                       [ y = 0            ]
    // 2) at or below and to the right of the line from (-3,0) to (0,4)  [ y = 4/3 * x + 4  ]
    // 3) at or below and to the left of the line from (3,0) to (0,4)    [ y = -4/3 * x + 4 ]
    if (y >= 0 && (y <= (4.0/3.0)*x + 4) && (y <= (-4.0/3.0)*x + 4))
    {
        return true;
    }

    return false;
}


// Defines two segments crossing
struct cross_res
{
    double x;
    double y;
    bool state;
};
// Code is based of https://ru.wikipedia.org/wiki/Пересечение_(евклидова_геометрия) and uses Cramer's rule
cross_res get_segment_crossing(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    // Determinant in Cramer's rule
    double det = (x4 - x3)*(y2 - y1) - (x2 - x1)*(y4 - y3);
    // System solution (yields s0, t0 we seek)
    double s0 = ((x4 - x3)*(y3 - y1) - (x3 - x1)*(y4 - y3)) / det;
    double t0 = ((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)) / det;

    // Check if crossing exists
    if ((0 < s0) && (0 < t0) && (s0 < 1) && (t0 < 1))
    {
        return cross_res{x1 + s0*(x2 - x1), y1 + s0*(y2 - y1), true};
    }
    else
    {
        return cross_res{0, 0, false};
    }
}


// Check if segment crosses D and get the crossing point if so
cross_res get_segment_crossing_with_D(double x1, double y1, double x2, double y2)
{
    cross_res crossing = get_segment_crossing(-3, 0, 0, 4, x1, y1, x2, y2);
    // Checks segment [(-3,0), (0, 4)]
    if (crossing.state)
    {
        return crossing;
    }
    // Checks segment [(0,4), (3, 0)]
    crossing = get_segment_crossing(0, 4, 3, 0, x1, y1, x2, y2);
    if (crossing.state)
    {
        return crossing;
    }
    // Checks segment [(-3,0), (3, 0)]
    if ((y1 == 0) && (y2 == 0)) // Segment might be at y = 0
    {
        if ((x1 < -3) || (x2 < -3))
        {
            crossing = {-3, 0, true};
            return crossing;
        }
        if ((x1 > 3) || (x2 > 3))
        {
            crossing = {3, 0, true};
            return crossing;
        }
    }
    crossing = get_segment_crossing(-3, 0, 3, 0, x1, y1, x2, y2);
    if (crossing.state)
    {
        return crossing;
    }

    // Otherwise there is no crossing
    crossing = {0, 0, false};
    return crossing;
}


// Differential operator over grid function
#define DIFF_OPER_A(func) \
    (-(a_ij[j*sizeM + i + 1]*(func[j*sizeM + i + 1] - func[j*sizeM + i])/h_x - a_ij[j*sizeM + i]*(func[j*sizeM + i] - func[j*sizeM + i - 1])/h_x)/h_x \
     -(b_ij[(j+1)*sizeM + i]*(func[(j+1)*sizeM + i] - func[j*sizeM + i])/h_y - b_ij[j*sizeM + i]*(func[j*sizeM + i] - func[(j-1)*sizeM + i])/h_y)/h_y)


// Returns value of the right part in point (x,y)
double f(double x, double y)
{
    return 1;
}


// Saves matrix to CSV file
bool save_matrix(string name, double *matrix, int N, int M)
{
    // Open file stream
    ofstream result_f;
    result_f.open("results/" + name + ".csv");
    // Check stream is open
    if (result_f.is_open())
    {
        // Fill in the CSV
        for (int j = 0; j < N; j++)
        {
            for (int i = 0; i < M; i++)
            {
                result_f << matrix[j*M + i];
                // Do not add ';' in the end
                if (i != M-1)
                {
                    result_f << ';';
                }
            }
            result_f << endl;
        }
        result_f.close();

        return true;
    }
    else
    {
        cout << "Unable to open file \"" << name + ".csv" << "\"!\n";
        return false;
    }
}


#define SYNC_DATA_BORDER(data) \
/* Sync direction up */ \
if (pid_j != 0) \
{ \
    MPI_Send(data + sizeM + 1, sizeM-2, MPI_DOUBLE, (pid_j-1)*pid_i_max + pid_i, 1, MPI_COMM_WORLD); \
} \
if (pid_j != pid_j_max-1) \
{ \
    MPI_Recv(data + (sizeN-1)*sizeM + 1, sizeM-2, MPI_DOUBLE, (pid_j+1)*pid_i_max + pid_i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
} \
/* Sync direction right */ \
if (pid_i != pid_i_max-1) \
{ \
    _Pragma("omp parallel for") \
    for (int j = 1; j < sizeN-1; j++) \
    { \
        column_buff[j-1] = data[j*sizeM + sizeM-2]; \
    } \
    MPI_Send(column_buff, sizeN-2, MPI_DOUBLE, (pid_j)*pid_i_max + pid_i + 1, 2, MPI_COMM_WORLD); \
} \
if (pid_i != 0) \
{ \
    MPI_Recv(column_buff, sizeN-2, MPI_DOUBLE, (pid_j)*pid_i_max + pid_i - 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
    _Pragma("omp parallel for") \
    for (int j = 1; j < sizeN-1; j++) \
    { \
        data[j*sizeM] = column_buff[j-1]; \
    } \
} \
/* Sync direction down */ \
if (pid_j != pid_j_max-1) \
{ \
    MPI_Send(data + (sizeN-2)*sizeM + 1, sizeM-2, MPI_DOUBLE, (pid_j+1)*pid_i_max + pid_i, 3, MPI_COMM_WORLD); \
} \
if (pid_j != 0) \
{ \
    MPI_Recv(data + 1, sizeM-2, MPI_DOUBLE, (pid_j-1)*pid_i_max + pid_i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
} \
/* Sync direction left */ \
if (pid_i != 0) \
{ \
    _Pragma("omp parallel for") \
    for (int j = 1; j < sizeN-1; j++) \
    { \
        column_buff[j-1] = data[j*sizeM + 1]; \
    } \
    MPI_Send(column_buff, sizeN-2, MPI_DOUBLE, (pid_j)*pid_i_max + pid_i - 1, 4, MPI_COMM_WORLD); \
} \
if (pid_i != pid_i_max-1) \
{ \
    MPI_Recv(column_buff, sizeN-2, MPI_DOUBLE, (pid_j)*pid_i_max + pid_i + 1, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE); \
    _Pragma("omp parallel for") \
    for (int j = 1; j < sizeN-1; j++) \
    { \
        data[j*sizeM + sizeM-1] = column_buff[j-1]; \
    } \
}


// #################################### //
int main(int argc, char *argv[])
{
    // Get basic MPI info
    int pid, total_pids;
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &total_pids);
    // Grid coordinates
    int dims[2] = {0,0};
    MPI_Dims_create(total_pids, 2, dims);
    int pid_j_max = dims[0],  pid_i_max = dims[1];
    int pid_j = pid/pid_i_max, pid_i = pid%pid_i_max;
    string pid_str = "(" + to_string(pid_j) + "," + to_string(pid_i) + ")";
    string pid_debug_str = "PID" + pid_str + ": ";

    // Debug
    printf("%sProgram is now alive!\n", pid_debug_str.c_str());
    chrono::steady_clock::time_point prog_stared = chrono::steady_clock::now();

    // Check args
    if (argc < 4)
    {
        printf("%sN,M, eps, delta are not set!", pid_debug_str.c_str());
        return -1;
    }

    string exec_name = argv[0];
    // Grid setup
    int N = stoi(argv[1]);
    int M = stoi(argv[2]);
    // Delta
    double delta = stod(argv[3]);

    // Debug
    printf("%sReceived parameters: %d, %d, %f\n", pid_debug_str.c_str(), N, M, delta);
    chrono::steady_clock::time_point solve_stared = chrono::steady_clock::now();

    // Max and min values of coordinates (define rectangle)
    double xmin = -3, xmax = 3, ymin = 0, ymax = 4;
    // Grid steps
    double h_x = (xmax - xmin) / M, h_y = (ymax - ymin) / N;

    // Epsilon = max of grid steps
    double eps = max(h_x, h_y);
    eps *= eps;

    // Debug
    printf("%sOther parameters init: Success\n", pid_debug_str.c_str());

    // Debug
    printf("%sMPI init: Success. (Pros %d of %d in total)\n", pid_debug_str.c_str(), pid, total_pids);

    // Matrix sizes
    int sizeN = (N/pid_j_max + 1 + (pid_j==(pid_j_max-1) ? N%pid_j_max : 0));
    int sizeM = (M/pid_i_max + 1 + (pid_i==(pid_i_max-1) ? M%pid_i_max : 0));
    int size = sizeN*sizeM;
    // Allocalte a_ij, b_ij and F_ij
    double *a_ij = reinterpret_cast<double*>(malloc(size * sizeof(double)));
    double *b_ij = reinterpret_cast<double*>(malloc(size * sizeof(double)));
    double *F_ij = reinterpret_cast<double*>(malloc(size * sizeof(double)));
    // Allocalte w_ij_prev, w_ij_curr
    double *w_ij_prev = reinterpret_cast<double*>(malloc(size * sizeof(double)));
    double *w_ij_curr = reinterpret_cast<double*>(malloc(size * sizeof(double)));
    // Allocate r_ij
    double *r_ij = reinterpret_cast<double*>(malloc(size * sizeof(double)));
    // Allocate column buffer for MPI
    double *column_buff = reinterpret_cast<double*>(malloc((sizeN-2) * sizeof(double)));

    // Debug
    printf("%sAllocating memory: Success\n", pid_debug_str.c_str());

    // Init w_ij_prev as zero matrix
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < sizeN; j++)
    {
        for (int i = 0; i < sizeM; i++)
        {
            w_ij_prev[j*sizeM + i] = 0;
            r_ij[j*sizeM + i] = 0;
        }
    }

    // Debug
    printf("%sInit w_ij: Success\n", pid_debug_str.c_str());

    // Fill in a_ij, b_ij and F_ij
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < sizeN; j++)
    {
        for (int i = 0; i < sizeM; i++)
        {
            // Half coordinates
            double x_i_mhalf = xmin + ((i + pid_i * M/pid_i_max - (pid_i?1:0)) + 0.5)*h_x; // i+1 shifts our triangle properly
            double x_i_phalf = xmin + ((i + pid_i * M/pid_i_max - (pid_i?1:0)) + 1.5)*h_x;
            double y_j_mhalf = ymax - ((j + pid_j * N/pid_j_max - (pid_j?1:0)) - 0.5)*h_y; // "ymax -" instead of "ymin + "
            double y_j_phalf = ymax - ((j + pid_j * N/pid_j_max - (pid_j?1:0)) + 0.5)*h_y; //  inverts triangle upwards

            // a_ij
            // Inside D
            if (is_point_in_D(x_i_mhalf, y_j_mhalf) && is_point_in_D(x_i_mhalf, y_j_phalf))
            {
                a_ij[j*sizeM + i] = 1;
            }
            // Outside D
            else if (!is_point_in_D(x_i_mhalf, y_j_mhalf) && !is_point_in_D(x_i_mhalf, y_j_phalf))
            {
                a_ij[j*sizeM + i] = 1/eps;
            }
            // Partially in D
            else
            {
                // Get point in D
                double point_in_D_y;
                if (is_point_in_D(x_i_mhalf, y_j_mhalf))
                {
                    point_in_D_y = y_j_mhalf;
                }
                else if (is_point_in_D(x_i_mhalf, y_j_phalf))
                {
                    point_in_D_y = y_j_phalf;
                }

                // Calculate length of the segment part in D
                cross_res crossing = get_segment_crossing_with_D(x_i_mhalf, y_j_mhalf, x_i_mhalf, y_j_phalf);
                double l_ij = sqrt(pow(crossing.x - x_i_mhalf, 2) + pow(crossing.y - point_in_D_y, 2));

                // According to formula
                a_ij[j*sizeM + i] = l_ij / h_y + (1 - l_ij / h_y) / eps;
            }

            // b_ij
            // Inside D
            if (is_point_in_D(x_i_mhalf, y_j_mhalf) && is_point_in_D(x_i_phalf, y_j_mhalf))
            {
                b_ij[j*sizeM + i] = 1;
            }
            // Outside D
            else if (!is_point_in_D(x_i_mhalf, y_j_mhalf) && !is_point_in_D(x_i_phalf, y_j_mhalf))
            {
                b_ij[j*sizeM + i] = 1/eps;
            }
            // Partially in D
            else
            {
                // Get point in D
                double point_in_D_x;
                if (is_point_in_D(x_i_mhalf, y_j_mhalf))
                {
                    point_in_D_x = x_i_mhalf;
                }
                else if (is_point_in_D(x_i_phalf, y_j_mhalf))
                {
                    point_in_D_x = x_i_phalf;
                }

                // Calculate length of the segment part in D
                cross_res crossing = get_segment_crossing_with_D(x_i_mhalf, y_j_mhalf, x_i_phalf, y_j_mhalf);
                double l_ij = sqrt(pow(crossing.x - point_in_D_x, 2) + pow(crossing.y - y_j_mhalf, 2));

                // According to formula
                b_ij[j*sizeM + i] = l_ij / h_x + (1 - l_ij / h_x) / eps;
            }

            // F_ij
            // Inside D
            if (is_point_in_D(x_i_mhalf, y_j_mhalf) && is_point_in_D(x_i_mhalf, y_j_phalf) &&
                is_point_in_D(x_i_phalf, y_j_mhalf) && is_point_in_D(x_i_phalf, y_j_phalf))
            {
                F_ij[j*sizeM + i] = f(xmin + i*h_x, ymin + j*h_y);
            }
            // Outside D
            else if (!is_point_in_D(x_i_mhalf, y_j_mhalf) && !is_point_in_D(x_i_mhalf, y_j_phalf) &&
                     !is_point_in_D(x_i_phalf, y_j_mhalf) && !is_point_in_D(x_i_phalf, y_j_phalf))
            {
                F_ij[j*sizeM + i] = 0;
            }
            // Partially in D
            else
            {
                // Points, that lie in П_ij {E D
                double points[7][2];
                int points_count = 0;

                // F(x, y) somewhere in П_ij {E D
                double f_in_S_ij;
                // Set f_in_S_ij and gather all points in П_ij {E D
                // (-,-)
                if (is_point_in_D(x_i_mhalf, y_j_mhalf))
                {
                    f_in_S_ij = f(x_i_mhalf, y_j_mhalf);
                    points[points_count][0] = x_i_mhalf;
                    points[points_count][1] = y_j_mhalf;
                    points_count++;
                }
                // (-,-) -> (+,-)
                cross_res crossing = get_segment_crossing_with_D(x_i_mhalf, y_j_mhalf, x_i_phalf, y_j_mhalf);
                if (crossing.state)
                {
                    points[points_count][0] = crossing.x;
                    points[points_count][1] = crossing.y;
                    points_count++;
                }
                // (+,-)
                if (is_point_in_D(x_i_phalf, y_j_mhalf))
                {
                    f_in_S_ij = f(x_i_phalf, y_j_mhalf);
                    points[points_count][0] = x_i_phalf;
                    points[points_count][1] = y_j_mhalf;
                    points_count++;
                }
                // (+,-) -> (+,+)
                crossing = get_segment_crossing_with_D(x_i_phalf, y_j_mhalf, x_i_phalf, y_j_phalf);
                if (crossing.state)
                {
                    points[points_count][0] = crossing.x;
                    points[points_count][1] = crossing.y;
                    points_count++;
                }
                // (+,+)
                if (is_point_in_D(x_i_phalf, y_j_phalf))
                {
                    f_in_S_ij = f(x_i_phalf, y_j_phalf);
                    points[points_count][0] = x_i_phalf;
                    points[points_count][1] = y_j_phalf;
                    points_count++;
                }
                // (+,+) -> (-,+)
                crossing = get_segment_crossing_with_D(x_i_phalf, y_j_phalf, x_i_mhalf, y_j_phalf);
                if (crossing.state)
                {
                    points[points_count][0] = crossing.x;
                    points[points_count][1] = crossing.y;
                    points_count++;
                }
                // (-,+)
                if (is_point_in_D(x_i_mhalf, y_j_phalf))
                {
                    f_in_S_ij = f(x_i_mhalf, y_j_phalf);
                    points[points_count][0] = x_i_mhalf;
                    points[points_count][1] = y_j_phalf;
                    points_count++;
                }
                // (-,+) -> (-,-)
                crossing = get_segment_crossing_with_D(x_i_mhalf, y_j_phalf, x_i_mhalf, y_j_mhalf);
                if (crossing.state)
                {
                    points[points_count][0] = crossing.x;
                    points[points_count][1] = crossing.y;
                    points_count++;
                }

                // Repeat the first point in the end
                points[points_count][0] = points[0][0];
                points[points_count][1] = points[0][1];

                // https://ru.wikihow.com/найти-площадь-многоугольника
                double a = 0;
                for (int k = 0; k < points_count; k++)
                {
                    a += points[k][1] * points[k + 1][0];
                }
                double b = 0;
                for (int k = 0; k < points_count; k++)
                {
                    b += points[k][0] * points[k + 1][1];
                }
                double S_ij = abs(a - b) / 2;
                
                // According to formula
                F_ij[j*sizeM + i] = S_ij * f_in_S_ij / (h_x * h_y);
            }
        }
    }

    // Debug
    printf("%sInit a_ij, b_ij, F_ij: Success\n", pid_debug_str.c_str());
    // Save matrices
    save_matrix(exec_name + "_p" + pid_str + "_a", a_ij, sizeN, sizeM);
    save_matrix(exec_name + "_p" + pid_str + "_b", b_ij, sizeN, sizeM);
    save_matrix(exec_name + "_p" + pid_str + "_f", F_ij, sizeN, sizeM);

    bool loop_cond = true;
    double iter_delta;
    // We must count iterations
    int iters = 0;
    // Iterations
    while (loop_cond)
    {
        // Compute r_ij
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < sizeN-1; j++)
        {
            for (int i = 1; i < sizeM-1; i++)
            {
                r_ij[j*sizeM + i] = DIFF_OPER_A(w_ij_prev) - F_ij[j*sizeM + i];
            }
        }
        SYNC_DATA_BORDER(r_ij)

        // Compute iteration parameter (division of dot products)
        double dot_product1 = 0;
        double dot_product2 = 0;
        #pragma omp parallel for reduction(+:dot_product1,dot_product2)
        for (int j = 1; j < sizeN-1; j++)
        {
            for (int i = 1; i < sizeM-1; i++)
            {
                dot_product1 += r_ij[j*sizeM + i] * r_ij[j*sizeM + i];
                dot_product2 += DIFF_OPER_A(r_ij) * r_ij[j*sizeM + i];
            }
        }
        MPI_Allreduce(&dot_product1, &dot_product1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&dot_product2, &dot_product2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double iter_param = dot_product1/dot_product2;

        // Compute w_ij_curr
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < sizeN-1; j++)
        {
            for (int i = 1; i < sizeM-1; i++)
            {
                w_ij_curr[j*sizeM + i] = w_ij_prev[j*sizeM + i] - iter_param*r_ij[j*sizeM + i];
            }
        }

        // Compute iteration delta (Euclidean norm)
        iter_delta = 0;
        double diff;
        #pragma omp parallel for private(diff) reduction(+:iter_delta)
        for (int j = 1; j < sizeN-1; j++)
        {
            for (int i = 1; i < sizeM-1; i++)
            {
                diff = w_ij_curr[j*sizeM + i] - w_ij_prev[j*sizeM + i];
                iter_delta += diff*diff;
            }
        }
        MPI_Allreduce(&iter_delta, &iter_delta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        iter_delta = sqrt(iter_delta);

        // Depending on stop condition
        if (iter_delta < delta)
        {
            // Stop loop
            loop_cond = false;
        }
        else
        {
            // Copy current w_ij to previos w_ij
            #pragma omp parallel for collapse(2)
            for (int j = 1; j < sizeN-1; j++)
            {
                for (int i = 1; i < sizeM-1; i++)
                {
                    w_ij_prev[j*sizeM + i] = w_ij_curr[j*sizeM + i];
                }
            }
            SYNC_DATA_BORDER(w_ij_prev)
        }

        iters++;
    }

    // Debug
    printf("%sLoop edned at iteration %d with delta %f\n", pid_debug_str.c_str(), iters, iter_delta);
    chrono::steady_clock::time_point solve_ended = chrono::steady_clock::now();

    // Save results
    save_matrix(exec_name + "_p" + pid_str + "_res", w_ij_curr, sizeN, sizeM);

    // Debug
    printf("%sSaving results: Success\n", pid_debug_str.c_str());

    // Free allocated memory
    free(a_ij);
    free(b_ij);
    free(F_ij);
    free(w_ij_prev);
    free(w_ij_curr);
    free(r_ij);
    free(column_buff);

    // Debug
    printf("%sFreeing memory: Success\n", pid_debug_str.c_str());

    // Debug
    chrono::steady_clock::time_point prog_ended = chrono::steady_clock::now();
    printf("%sSolver run time: %ld.%ld[s]\n", pid_debug_str.c_str(),
           chrono::duration_cast<std::chrono::milliseconds>(solve_ended - solve_stared).count()/1000,
           chrono::duration_cast<std::chrono::milliseconds>(solve_ended - solve_stared).count()%1000);
    printf("%sProgramm run time: %ld.%ld[s]\n", pid_debug_str.c_str(),
           chrono::duration_cast<std::chrono::milliseconds>(prog_ended - prog_stared).count()/1000,
           chrono::duration_cast<std::chrono::milliseconds>(prog_ended - prog_stared).count()%1000);

    // MPI ending
    MPI_Finalize();

    return 0;
}