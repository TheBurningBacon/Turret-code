
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

using namespace std;

// returns {time taken,vx when hit, vy when hit}
// time taken is -1 if it missed
vector<float> calc(float v_i, float angle,float target_x, float target_y,float dt=0.1) {
    float g=9.81; //gravity
    float k=1; //drag constant (PcA/2m)
    
    float x=0;
    float y=0;
    float vx = v_i * cos(angle);
    float vy = v_i * sin(angle);
    float t = 0;

    float time_to_target=-1;
    float tolerance = 0.05; //Distance tolerance to condier target hit

    vector<float> time_history;
    vector<float> x_v_history;
    vector<float> y_v_history;
    vector<float> x_pos_history;
    vector<float> y_pos_history;

    while ((y>target_y || vy>0) && x>target_x) {
        float v = sqrt((vx*vx) + (vy*vy));

        //float drag_magnitude = k * v^2

        //float ax_drag = -k * v^2 * (vx / v) = -k * v * vx
        //float ay_drag = -k * v^2 * (vy / v) = -k * v * vy
        
        float ax = -k * v *vx;
        float ay = -g - (k * v *vy);

        vx = vx + (ax * dt);
        vy = vy + (ay * dt);

        x = x + (vx * dt);
        y = y + (vy * dt);

        t = t + dt;

        time_history.push_back(t);
        x_v_history.push_back(vx);
        y_v_history.push_back(vy);
        x_pos_history.push_back(x);
        y_pos_history.push_back(y);

        if (time_to_target==-1 &&
            abs(x - target_x) <= tolerance &&
            abs(y - target_y) <= tolerance) 
            {
            time_to_target = t;
            break;
            }
    }

    return {time_to_target,vx,vy};
}

struct TurretAngles {
    bool found = false;
    float horizontal_angle = 0;
    float vertical_angle = 0;
    float horizontal_velocity = 0;
};

bool height_error_at_target(float initial_velocity, float pitch, float target_x,
                            float target_y, float dt, float& height_error) {
    constexpr float g = 9.81f;
    constexpr float k = 1.0f;
    float x = 0;
    float y = 0;
    float vx = initial_velocity * cos(pitch);
    float vy = initial_velocity * sin(pitch);

    for (int step = 0; step < 200000 && vx > 0; ++step) {
        float previous_x = x;
        float previous_y = y;
        float speed = sqrt(vx * vx + vy * vy);
        float ax = -k * speed * vx;
        float ay = -g - k * speed * vy;

        vx += ax * dt;
        vy += ay * dt;
        x += vx * dt;
        y += vy * dt;

        if (x >= target_x) {
            float fraction = (target_x - previous_x) / (x - previous_x);
            float y_at_target = previous_y + fraction * (y - previous_y);
            height_error = y_at_target - target_y;
            return true;
        }
    }

    return false;
}

TurretAngles find_turret_angles(float distance, float horizontal_angle,
                               float target_vertical_angle, float initial_velocity,
                               float dt = 0.01f) {
    TurretAngles result;
    constexpr float pi = 3.14159265358979323846f;

    if (distance <= 0 || initial_velocity <= 0 || dt <= 0 ||
        !isfinite(distance) || !isfinite(horizontal_angle) ||
        !isfinite(target_vertical_angle) || !isfinite(initial_velocity) ||
        !isfinite(dt) || abs(target_vertical_angle) >= pi / 2) {
        return result;
    }

    float target_x = distance * cos(target_vertical_angle);
    float target_y = distance * sin(target_vertical_angle);
    float pitch_min = -pi / 2 + 0.01f;
    float pitch_max = pi / 2 - 0.01f;
    constexpr int samples = 720;
    float previous_pitch = pitch_min;
    float previous_error = 0;
    bool previous_valid = false;

    for (int sample = 0; sample <= samples; ++sample) {
        float pitch = pitch_min + (pitch_max - pitch_min) * sample / samples;
        float error = 0;
        bool valid = height_error_at_target(initial_velocity, pitch, target_x,
                                            target_y, dt, error);

        if (valid && abs(error) <= 0.0001f) {
            result.found = true;
            result.horizontal_angle = horizontal_angle;
            result.vertical_angle = pitch;
            result.horizontal_velocity = initial_velocity * cos(pitch);
            return result;
        }

        if (valid && previous_valid &&
            ((error < 0 && previous_error > 0) ||
             (error > 0 && previous_error < 0))) {
            float low = previous_pitch;
            float high = pitch;
            float low_error = previous_error;

            for (int iteration = 0; iteration < 35; ++iteration) {
                float middle = (low + high) / 2;
                float middle_error = 0;
                if (!height_error_at_target(initial_velocity, middle, target_x,
                                            target_y, dt, middle_error)) {
                    break;
                }

                if ((middle_error < 0) == (low_error < 0)) {
                    low = middle;
                    low_error = middle_error;
                } else {
                    high = middle;
                }
            }

            float pitch_solution = (low + high) / 2;
            result.found = true;
            result.horizontal_angle = horizontal_angle;
            result.vertical_angle = pitch_solution;
            result.horizontal_velocity = initial_velocity * cos(pitch_solution);
            return result;
        }

        previous_pitch = pitch;
        previous_error = error;
        previous_valid = valid;
    }

    return result;
}

class Target {
private:
    string id;
    //pos,velocity,acceleration
    float x;
    float y;
    float z;
    vector<tuple<float,float,float>> pos_history;
    float vx=0;
    float vy=0;
    float vz=0;
    vector<tuple<float,float,float>> vel_history;
    float ax=0;
    float ay=0;
    float az=0;
    
public:
    Target(float x,float y,float z,string id): x{x}, y{y}, z{z}, id{id}{
        
    }

    string get_id() {
        return id;}

    tuple<float,float,float> cacl_pos_at_time(float time) {
        float kv = 0;//These need to be calibrated
        float ka = 0;

        float new_x=x+(vx*time*kv)+(ax*time*ka);
        float new_y=y+(vy*time*kv)+(ay*time*ka);
        float new_z=z+(vz*time*kv)+(az*time*ka);

        return {new_x,new_y,new_z};

    }

    void push_data(float new_x,float new_y,float new_z) {

        int frames_checked=10;

        pos_history.push_back({new_x,new_y,new_z});


        //gets avarage change from last ten frames

        float old_old_x;
        float old_old_y;
        float old_old_z;
        float avg_v_x=0;
        float avg_v_y=0;
        float avg_v_z=0;
        for (int i=0;i<min((int)pos_history.size(),frames_checked);i++) {
            //old is the value of the current loop
            //old_old is the value of the previous loop
            float old_x=get<0>(pos_history[i]);
            float old_y=get<1>(pos_history[i]);
            float old_z=get<2>(pos_history[i]);

            //skips the first loop and defines old_old's
            if (i!=0) {
                avg_v_x+=old_x-old_old_x;
                avg_v_y+=old_y-old_old_y;
                avg_v_z+=old_z-old_old_z;

            }
            else {
                old_old_x=old_x;
                old_old_y=old_y;
                old_old_z=old_z;
            }
        }

        avg_v_x/=min((int)pos_history.size(),frames_checked);
        avg_v_y/=min((int)pos_history.size(),frames_checked);
        avg_v_z/=min((int)pos_history.size(),frames_checked);

        vel_history.push_back({avg_v_x,avg_v_y,avg_v_z});


        float old_old_vx;
        float old_old_vy;
        float old_old_vz;
        float avg_a_x=0;
        float avg_a_y=0;
        float avg_a_z=0;
        for (int i=0;i<min((int)vel_history.size(),frames_checked);i++) {
            //old is the value of the current loop
            //old_old is the value of the previous loop
            float old_vx=get<0>(vel_history[i]);
            float old_vy=get<1>(vel_history[i]);
            float old_vz=get<2>(vel_history[i]);

            //skips the first loop and defines old_old's
            if (i!=0) {
                avg_a_x+=old_vx-old_old_vx;
                avg_a_y+=old_vy-old_old_vy;
                avg_a_z+=old_vz-old_old_vz;

            }
            else {
                old_old_vx=old_vx;
                old_old_vy=old_vy;
                old_old_vz=old_vz;
            }
        }

        avg_v_x/=min((int)vel_history.size(),frames_checked);
        avg_v_y/=min((int)vel_history.size(),frames_checked);
        avg_v_z/=min((int)vel_history.size(),frames_checked);





        x=new_x;
        y=new_y;
        z=new_z;
        vx=avg_v_x;
        vy=avg_v_y;
        vz=avg_v_z;
        ax=avg_a_x;
        ay=avg_a_y;
        az=avg_a_z;
    }

};


int main() {

    TurretAngles angles = find_turret_angles(0.1f, 0.5f, 0.0f, 5.0f);
    if (angles.found) {
        cout << "Yaw: " << angles.horizontal_angle
             << ", pitch: " << angles.vertical_angle
             << ", horizontal velocity: " << angles.horizontal_velocity << '\n';
    } else {
        cout << "Target is out of range\n";
    }

    return 0;
}