// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }


//Making the structure who will be given as parameter
typedef struct {

    int nrthread;
    int id;    
    pthread_barrier_t *barrier;    
    ppm_image *image;
    ppm_image *new_image;  
    char argim;
    int step_x;
    int step_y;
    unsigned char sigma;
    unsigned char **grid;
    ppm_image **contour_map;

} threaduri;

int min(int x, int y){
    return(x < y) ? x : y;
}

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "../checker/contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);
}


void *function(void *arg) {

//Rescale function parallelization

    threaduri *th = (threaduri*)arg;
    int thread_id = th->id;

    //Calculate the parts of the image who will be rescaled by the given thread

    int start = thread_id * (double)th->new_image->x / th->nrthread;
    int end = min((thread_id + 1) * (double)th->new_image->x / th->nrthread, th->new_image->x);
    uint8_t sample[3];

    // we only rescale downwards
    if (th->image->x > RESCALE_X || th->image->y > RESCALE_Y) {
        // Use bicubic interpolation for scaling
        for (int i = start; i < end; i++) {
            for (int j = 0; j < th->new_image->y; j++) {
                float u = (float)i / (float)(th->new_image->x - 1);
                float v = (float)j / (float)(th->new_image->y - 1);
                sample_bicubic(th->image, u, v, sample);

                th->new_image->data[i * th->new_image->y + j].red = sample[0];
                th->new_image->data[i * th->new_image->y + j].green = sample[1];
                th->new_image->data[i * th->new_image->y + j].blue = sample[2];
            }
        }
    } else {
        th->new_image = th->image;
    }

    //waiting for all the threads to rescale the image
    pthread_barrier_wait(th->barrier);

    //sample_grid function parallelization

    int p = th->new_image->x / th->step_x;
    int q = th->new_image->y / th->step_y;

    start = thread_id * (double)p / th->nrthread;
    end = min((thread_id + 1) * (double)p / th->nrthread, p);


    for (int i = start; i < end; i++) {
        for (int j = 0; j < q; j++) {
            ppm_pixel curr_pixel = th->new_image->data[i * th->step_x * th->new_image->y + j * th->step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > th->sigma) {
                th->grid[i][j] = 0;
            } else {
                th->grid[i][j] = 1;
            }
        }
    }

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them

    // Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
    // Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
    // pixel values compare to the `sigma` reference value. The points are taken at equal distances
    // in the original image, based on the `step_x` and `step_y` arguments.

    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = th->new_image->data[i * th->step_x * th->new_image->y + th->new_image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > th->sigma) {
            th->grid[i][q] = 0;
        } else {
            th->grid[i][q] = 1;
        }
    }
    for (int j = start; j < end; j++) {
        ppm_pixel curr_pixel = th->new_image->data[(th->new_image->x - 1) * th->new_image->y + j * th->step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > th->sigma) {
            th->grid[p][j] = 0;
        } else {
            th->grid[p][j] = 1;
        }
    }

    //parallelization of march function
    // Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
    // type of contour which corresponds to each subgrid. It determines the binary value of each
    // sample fragment of the original image and replaces the pixels in the original image with
    // the pixels of the corresponding contour image accordingly.

    for (int i = start; i < end; i++) 
    {
        for (int j = 0; j < q; j++) {
                unsigned char k = 8 * th->grid[i][j] + 4 * th->grid[i][j + 1] + 2 * th->grid[i + 1][j + 1] + 1 * th->grid[i + 1][j];
                update_image(th->new_image, th->contour_map[k], i * th->step_x, j * th->step_y);
	 }
    }

    pthread_exit(NULL);
}


int main(int argc, char *argv[]) {

    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    // 0. Initialize contour map
    ppm_image **contour_map = init_contour_map();

    //initialize the elements of the structure
    //the parameters of the parallelized functions

    ppm_image *image = read_ppm(argv[1]);
    int step_x = STEP;
    int step_y = STEP;
    int i, r;
    int nr = atoi(argv[3]);
    
    // making the vector of threads
    threaduri th[nr];

    //alloc memory for the new image
    th[0].new_image = (ppm_image *)malloc(sizeof(ppm_image));
    th[0].new_image->x = RESCALE_X;
    th[0].new_image->y = RESCALE_Y;

    //alloc memory for the data of the new image
    th[0].new_image->data = (ppm_pixel*)malloc(RESCALE_X * RESCALE_Y * sizeof(ppm_pixel));

    //creating and initializing the barrier
    pthread_barrier_t barrier;

    pthread_barrier_init(&barrier, NULL, nr);
    pthread_t tid[nr];

    int p = th[0].new_image->x / step_x;
    int q = th[0].new_image->y / step_y;

    //alloc memory for the grid
    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    //initializing the structure elements
    for (i = 0; i < nr; i++) {
       th[i].barrier = &barrier;
        th[i].new_image = th[0].new_image;
        th[i].new_image->x = RESCALE_X;
        th[i].new_image->y = RESCALE_Y;
        th[i].new_image->data = th[0].new_image->data;
        th[i].nrthread = nr;
        th[i].id = i;
        th[i].image = image;
        th[i].step_x = step_x;
        th[i].step_y = step_y;
        th[i].sigma = SIGMA;
        th[i].grid = grid;
        th[i].contour_map = contour_map;
    }
    
    //creating the threads
	for (i = 0; i < nr; i++) {
		r = pthread_create(&tid[i], NULL, function, &th[i]);
		if (r) {
			printf("Eroare la crearea thread-ului %d\n", i);
			exit(-1);
		}
	}

    //joining the threads
	for (i = 0; i < nr; i++) {
		r = pthread_join(tid[i], NULL);

		if (r) {
			printf("Eroare la asteptarea thread-ului %d\n", i);
			exit(-1);
		}
	}

    // 1. Rescale the image
    ppm_image *scaled_image = th[nr - 1].new_image;

    // 4. Write output
    write_ppm(scaled_image, argv[2]);

    free_resources(scaled_image, th->contour_map, th->grid, step_x);

    pthread_barrier_destroy(&barrier);

    return 0;
}
