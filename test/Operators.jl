@testset "Operators" begin
    @testset "Fraunhofer" begin
        arr = randn(ComplexF32, (7, 6, 3)) 
        o2d, d2o = Fraunhofer(arr)
        ss = sqrt(size(arr, 1) * size(arr, 2))
        @test o2d(copy(arr)) ≈ fft(arr, (1,2)) ./ ss 
        @test d2o(copy(arr)) ≈ ifft(arr, (1,2)) .* ss 

    
         
        o2d, d2o = Fraunhofer(arr, fftshiftFlag=true)
        @test o2d(copy(arr)) ≈ fftshift(fft(arr, (1,2)), (1,2)) ./ ss 
        @test d2o(copy(arr)) ≈ ifft(ifftshift(arr, (1,2)), (1,2)) .* ss 
        
    end
end
