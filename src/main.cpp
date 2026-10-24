#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/SparseExtra>

#include <iostream>
#include <fstream>
#include <cstdint>
#include <algorithm>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace Eigen;

template <typename M>
void truncatedSVD(const BDCSVD<M>& svd, int k, M& C, M& D) {
    C = svd.matrixU().leftCols(k);
    D = svd.matrixV().leftCols(k);

    for(int i = 0; i < k; i++) {
        D.col(i) *= svd.singularValues()(i);
    }
}

bool saveMatrixToImage(const MatrixXd mat, const std::string& path) {
    Matrix<uint8_t, Dynamic, Dynamic, RowMajor> save_image = (mat*255).cwiseMax(0.0).cwiseMin(255.0).cast<uint8_t>();

  	return stbi_write_png(path.c_str(), mat.cols(), mat.rows(), 1, save_image.data(), mat.cols());
}

std::vector<double> nBiggestEigenvaluesSymmetricMatrix(MatrixXd m, unsigned n) {
    SelfAdjointEigenSolver<MatrixXd> eigensolver(m);

    if (eigensolver.info() != Eigen::Success) {
        std::cerr << "Error: Could not compute eigensolver!" << std::endl;
        abort();
    }

    std::vector<double> nBiggestEigenvalues;

    auto eigenvalues = eigensolver.eigenvalues();
    for(int i = 0; i < n; i++) {
        nBiggestEigenvalues.push_back(eigenvalues[m.cols() - 1 - i]);
    }

    return nBiggestEigenvalues;
}

int main(int argc, char *argv[]) {
    const char *input_image_path = argc > 1 ? argv[1] : "res/deer.jpg";

    // TASK 1

    // Load the image using stb_image
    int width, height, channels;
    uint8_t *image_data = stbi_load(input_image_path, &width, &height, &channels, 1); 
    if (!image_data){
        std::cerr << "Error: Could not load image " << input_image_path << std::endl;
        return 1;
    }

    // Prepare Eigen matrices for each image
    MatrixXd original_matrix(height, width);

    for (int i = 0; i < height; i++){
        for (int j = 0; j < width; j++){
            int index = (i * width + j);

            original_matrix(i, j) = static_cast<double>(image_data[index])/255.0;
        }
    }

    MatrixXd matrix_ata = original_matrix.transpose()*original_matrix;

    std::cout << "Task 1: ||A^T*A|| = " << matrix_ata.norm() << std::endl;

    // TASK 2
    auto nBiggestEigenvalues = nBiggestEigenvaluesSymmetricMatrix(matrix_ata, 2);
    
    std::cout << "Task 2: The 2 biggest eigenvalues of A are: " << sqrt(nBiggestEigenvalues[0]) << ", " << sqrt(nBiggestEigenvalues[1]) << std::endl;

    // TASK 3
    saveMarket(matrix_ata, "res/ATA.mtx");

    std::string command = "./eigen1 ATA.mtx eigvec.txt hist.txt -e pi -etol 1.0e-10";

	std::cout << std::endl << command << std::endl;
	int result = system(("cd res && " + command).c_str());

	if (result) {
		std::cout << "Failed command execution" << std::endl;
		return 1;
	}

    // TASK 4
    command = "./eigen1 ATA.mtx eigvec.txt hist.txt -e pi -etol 1.0e-10 -shift 867.85";

	std::cout << command << std::endl;
	result = system(("cd res && " + command).c_str());

	if (result) {
		std::cout << "Failed command execution" << std::endl;
		return 1;
	}

    // TASK 5
    BDCSVD<MatrixXd> svd(original_matrix, ComputeThinU | ComputeThinV);

    MatrixXd sigma = svd.singularValues().asDiagonal();
    
    std::cout << "Task 5: ||Σ|| = " << sigma.norm() << std::endl;

    // TASK 6
    MatrixXd C30, D30;
    MatrixXd C90, D90;

    truncatedSVD<MatrixXd>(svd, 30, C30, D30);
    truncatedSVD<MatrixXd>(svd, 90, C90, D90);

    std::cout << "Task 6: C30 nnz: " << C30.nonZeros() << std::endl;
    std::cout << "Task 6: D30 nnz: " << D30.nonZeros() << std::endl;
    std::cout << "Task 6: C90 nnz: " << C90.nonZeros() << std::endl;
    std::cout << "Task 6: D90 nnz: " << D90.nonZeros() << std::endl;

    // TASK 7
    MatrixXd A30 = C30*D30.transpose();
    std::string output_path = "res/A30.png";
    saveMatrixToImage(A30, output_path);
	std::cout << "Task 7: Image saved to " << output_path << std::endl;
    
    MatrixXd A90 = C90*D90.transpose();
    output_path = "res/A90.png";
    saveMatrixToImage(A90, output_path);
	std::cout << "Task 7: Image saved to " << output_path << std::endl;

    // TASK 8
    int n = 400;
    MatrixXd checkerboard(n,n);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if((i/40 + j/40) % 2 == 0) {
                checkerboard(i,j) = 0.0;
            }
            else {
                checkerboard(i,j) = 1.0;
            }
        }
    }

    std::cout << "Task 8: Norm = "
              << checkerboard.norm() << std::endl;

    // TASK 9
    MatrixXd noise = (60.0/255.0) * MatrixXd::Random(n,n);
    MatrixXd noisy = checkerboard + noise;

    output_path = "res/checknoisy.png";
    saveMatrixToImage(noisy, output_path);
	std::cout << "Task 9: Image saved to " << output_path << std::endl;

    // TASK 10
    Eigen::BDCSVD<Eigen::MatrixXd> svdNoisy(
        noisy, Eigen::ComputeThinU | Eigen::ComputeThinV);

    VectorXd sigmaNoisy = svdNoisy.singularValues();

    std::cout << "Task 10: First singular value = "
              << sigmaNoisy(0) << std::endl;

    std::cout << "Task 10: Second singular value = "
              << sigmaNoisy(1) << std::endl;

    // TASK 11
    MatrixXd C6 = svdNoisy.matrixU().leftCols(6);
    MatrixXd C12 = svdNoisy.matrixU().leftCols(12);

    MatrixXd D6 = svdNoisy.matrixV().leftCols(6);
    MatrixXd D12 = svdNoisy.matrixV().leftCols(12);

    for(int i = 0; i < 6; i++) {
        D6.col(i) *= sigmaNoisy(i);
    }

    for(int i = 0; i < 12; i++) {
        D12.col(i) *= sigmaNoisy(i);
    }

    std::cout << "Task 11: C6 size = "
              << C6.rows() << " x " << C6.cols() << std::endl;

    std::cout << "Task 11: D6 size = "
              << D6.rows() << " x " << D6.cols() << std::endl;

    std::cout << "Task 11: C12 size = "
              << C12.rows() << " x " << C12.cols() << std::endl;

    std::cout << "Task 11: D12 size = "
              << D12.rows() << " x " << D12.cols() << std::endl;

    // TASK 12
    MatrixXd A6 = C6 * D6.transpose();
    MatrixXd A12 = C12 * D12.transpose();

    output_path = "res/check6.png";
    saveMatrixToImage(A6, output_path);
	std::cout << "Task 12: Image saved to " << output_path << std::endl;
    
    output_path = "res/check12.png";
    saveMatrixToImage(A12, output_path);
	std::cout << "Task 12: Image saved to " << output_path << std::endl;

    return 0;
}