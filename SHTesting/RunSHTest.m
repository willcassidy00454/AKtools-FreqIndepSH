% Tests the spherical harmonics of AKtools by simulating third order IRs
% from a source to a rotating receiver
function RunSHTest()
    close all
    clear all

    sample_rate = 48000;
    bit_depth = 24;
    
    params_dir = "Simulation Parameters/";
    
    GenerateRIRs(params_dir, ...
        "Outputs/", ...
        sample_rate, ...
        bit_depth);

    figure
    for plot_idx = 1:9
        PlotSHDirectivity("Outputs/", 1, 17, plot_idx, 400, "Azimuth");
    end
    figure
    for plot_idx = 1:9
        PlotSHDirectivity("Outputs/", 18, 17, plot_idx, 400, "Elevation");
    end
end

function PlotSHDirectivity(read_dir, start_index, num_indices, channel_to_plot, trunc_length_samples, label)
    irs = zeros(num_indices,trunc_length_samples);
    maxima = zeros(num_indices,1);
    thetas = zeros(num_indices,1);

    for plot_index = 1:num_indices
        file_index = start_index + plot_index - 1;
        [ir, ~] = audioread(read_dir + "Rec Rotation "+ file_index +"_R1_S1.wav");
        irs(plot_index,:) = ir(1:trunc_length_samples, channel_to_plot);

        maxima(plot_index) = max(abs(irs(plot_index,:)),[],"all");
        thetas(plot_index) = (plot_index - 1) * 2 * pi / (num_indices - 1);
    end

    nexttile;
    polarplot(thetas,maxima);
    rlim([0 1.0])
    title(label + " - Spherical Harmonic " + (channel_to_plot - 1))
end

function GenerateRIRs(params_dir, output_dir, sample_rate, bit_depth)
    room_dims = readmatrix(params_dir + "room_dimensions.dat");
    alphas = readmatrix(params_dir + "absorption_coeffs.dat");
    
    src_coords = readmatrix(params_dir + "src_coords.dat");
    rec_coords = readmatrix(params_dir + "rec_coords.dat");
    
    src_rotations = readmatrix(params_dir + "src_rotations.dat");
    rec_rotations = readmatrix(params_dir + "rec_rotations.dat");
    
    src_directivities = string(readcell(params_dir + "src_directivities.csv"));
    rec_directivities = string(readcell(params_dir + "rec_directivities_3rd_order.csv"));

    for rec_rot_index = 1:size(rec_rotations, 1)
        for third_order_sh_output_index = 1:16
            ir = GenerateSrcToRecIRs(src_coords, ...
                rec_coords, ...
                src_rotations, ...
                rec_rotations(rec_rot_index, :), ...
                src_directivities, ...
                rec_directivities, ...
                room_dims, ...
                alphas, ...
                sample_rate, ...
                "Rec Rotation " + rec_rot_index, ...
                true, ...
                third_order_sh_output_index);

            third_order_irs(third_order_sh_output_index, rec_rot_index, 1:length(ir)) = ir;
        end
    end
    
    % normalise all
    third_order_irs = third_order_irs / max(abs(third_order_irs),[],'all');

    for rec_rot_index = 1:size(rec_rotations, 1)
        SaveIRsMultichannelReceivers(third_order_irs(:,rec_rot_index,:),sample_rate,bit_depth,output_dir,"Rec Rotation " + rec_rot_index);
    end
end