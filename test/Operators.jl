@testset "Operators" begin
    @testset "Fraunhofer" begin
        arr = randn(ComplexF32, (7, 6, 3)) 
        object2detector, detector2object = Fraunhofer(arr)
        ss = sqrt(size(arr, 1) * size(arr, 2))
        @test object2detector(copy(arr)) ≈ fft(arr, (1,2)) ./ ss 
        @test detector2object(copy(arr)) ≈ ifft(arr, (1,2)) .* ss 

    
         
        object2detector, detector2object = Fraunhofer(arr, fftshiftSwitch=true)
        @test object2detector(copy(arr)) ≈ fftshift(fft(arr, (1,2)), (1,2)) ./ ss 
        @test detector2object(copy(arr)) ≈ ifft(ifftshift(arr, (1,2)), (1,2)) .* ss 
        
    end
end
