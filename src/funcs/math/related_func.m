% 注意：MATLAB求得的旋转矩阵是我们平常习惯的转置！

% Rodrigues公式做旋转向量和旋转矩阵的转换
R = rotationVectorToMatrix(vec);
vec = rotationMatrixToVector(R);

% 四元数、旋转矩阵
rotm = quat2rotm(quat);
quat = rotm2quat(rotm);

% 
