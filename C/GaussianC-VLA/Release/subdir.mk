################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables
C_SRCS += \
../Timer.c \
../ThreadPool.c \
../gaussian.c \
../solve_system.c

C_DEPS += \
./Timer.d \
../ThreadPool.d \
./gaussian.d \
./solve_system.d

OBJS += \
./Timer.o \
../ThreadPool.o \
./gaussian.o \
./solve_system.o


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C Compiler'
	gcc -std=gnu11 -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean--2e-

clean--2e-:
	-$(RM) ./Timer.d ./Timer.o ./gaussian.d ./gaussian.o ./solve_system.d ./solve_system.o

.PHONY: clean--2e-

