/*
 * Software License Agreement (BSD License)
 *
 *  Point Cloud Library (PCL) - www.pointclouds.org
 *  Copyright (c) 2010, Willow Garage, Inc.
 *  Copyright (c) 2012-, Open Perception, Inc.
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id$
 *
 */
#ifndef PCL_REGISTRATION_IMPL_ICP_ANDERSON_HPP_
#define PCL_REGISTRATION_IMPL_ICP_ANDERSON_HPP_

#include <pcl/registration/boost.h>
#include <pcl/registration/exceptions.h>


///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar> typename pcl::AndersonIterativeClosestPoint<PointSource, PointTarget>::Vector6
pcl::AndersonIterativeClosestPoint<PointSource, PointTarget, Scalar>::Matrix42Vector6 (const Matrix4 m)
{
  Vector6 v;
  Matrix3 s = m.block(0,0,3,3);
  v.head(3) = s.eulerAngles(0, 1, 2);
  v.tail(3) = m.col(3).head(3);
  return v;
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar> typename pcl::AndersonIterativeClosestPoint<PointSource, PointTarget>::Matrix4
pcl::AndersonIterativeClosestPoint<PointSource, PointTarget, Scalar>::Vector62Matrix4 (const Vector6 v)
{
  Matrix3 s (Eigen::AngleAxis<Scalar>(v(0), Vector3::UnitX())
    * Eigen::AngleAxis<Scalar>(v(1), Vector3::UnitY())
    * Eigen::AngleAxis<Scalar>(v(2), Vector3::UnitZ()));
  Matrix4 m = Matrix4::Zero();
  m.block(0,0,3,3) = s;
  m(3,3) = 1;
  m.col(3).head(3) = v.tail(3);
  return m;
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar> int
pcl::AndersonIterativeClosestPoint<PointSource, PointTarget, Scalar>::alphas_cond (VectorX alphas)
{
  return alpha_limit_min_ < alphas.minCoeff() && alphas.maxCoeff() < alpha_limit_max_ && alphas(alphas.size()-1) > 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar> typename pcl::AndersonIterativeClosestPoint<PointSource, PointTarget>::VectorX
pcl::AndersonIterativeClosestPoint<PointSource, PointTarget, Scalar>::get_alphas_lstsq (const Matrix6X f)
{
  Matrix6X A = f.leftCols(f.cols()-1);
  A *= -1;
  A += f.rightCols(1) * VectorX::Constant(f.cols()-1, 1).transpose();
  VectorX sol = A.colPivHouseholderQr().solve(f.rightCols(1));
  sol.conservativeResize(sol.size()+1);
  sol[sol.size()-1] = 0;
  sol[sol.size()-1] = 1-sol.sum();
  return sol;
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar> typename pcl::AndersonIterativeClosestPoint<PointSource, PointTarget>::Vector6
pcl::AndersonIterativeClosestPoint<PointSource, PointTarget, Scalar>::get_next_u (const Matrix6X u, const Matrix6X g, const Matrix6X f)
{
  int i = 1;
  Vector6 sol = ((1-beta_)*u.col(u.cols()-1) + beta_*g.col(g.cols()-1));
  VectorX sol_alphas(1);
  sol_alphas << 1;

  for (i = 2; i <= f.cols(); i++)
  {
    VectorX alphas = get_alphas_lstsq(f.rightCols(i));
    if (!alphas_cond(alphas))
      break;
    sol = (1-beta_)*u.rightCols(i)*alphas + beta_*g.rightCols(i)*alphas;
    sol_alphas = alphas;
  }
  return sol;
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar> inline double
pcl::AndersonIterativeClosestPoint<PointSource, PointTarget, Scalar>::getMSE (Matrix4 fitness_transformation_)
{
  Matrix4 old = final_transformation_;
  final_transformation_ = fitness_transformation_;
  double res = convergence_criteria_->calculateMSE (*correspondences_);
  final_transformation_ = old;
  return res;
}

///////////////////////////////////////////////////////////////////////////////////////////
template <typename PointSource, typename PointTarget, typename Scalar> inline void
pcl::AndersonIterativeClosestPoint<PointSource, PointTarget, Scalar>::computeTransformation (PointCloudSource &output, const Eigen::Matrix4f& guess)
{
  double err = std::numeric_limits<double>::max ();
  double prev_err = std::numeric_limits<double>::max ();

  Matrix6X u(6,0), g(6,0), f(6,0);
  Vector6 u_next, u_k;

  // Point cloud containing the correspondences of each point in <input, indices>
  PointCloudSourcePtr input_transformed (new PointCloudSource);

  nr_iterations_ = 0;
  converged_ = false;
  bool snd_converged = false;

  // Initialise final transformation to the guessed one
  final_transformation_ = guess;

  // If the guessed transformation is non identity
  if (guess != Matrix4::Identity ())
  {
    input_transformed->resize (input_->size ());
     // Apply guessed transformation prior to search for neighbours
    transformCloud (*input_, *input_transformed, guess);
  }
  else
    *input_transformed = *input_;
 
  transformation_ = Matrix4::Identity ();

  // Make blobs if necessary
  determineRequiredBlobData ();
  PCLPointCloud2::Ptr target_blob (new PCLPointCloud2);
  if (need_target_blob_)
    pcl::toPCLPointCloud2 (*target_, *target_blob);

  // Pass in the default target for the Correspondence Estimation/Rejection code
  correspondence_estimation_->setInputTarget (target_);
  if (correspondence_estimation_->requiresTargetNormals ())
    correspondence_estimation_->setTargetNormals (target_blob);
  // Correspondence Rejectors need a binary blob
  for (size_t i = 0; i < correspondence_rejectors_.size (); ++i)
  {
    registration::CorrespondenceRejector::Ptr& rej = correspondence_rejectors_[i];
    if (rej->requiresTargetPoints ())
      rej->setTargetPoints (target_blob);
    if (rej->requiresTargetNormals () && target_has_normals_)
      rej->setTargetNormals (target_blob);
  }

  convergence_criteria_->setMaximumIterations (max_iterations_);
  convergence_criteria_->setRelativeMSE (euclidean_fitness_epsilon_);
  convergence_criteria_->setTranslationThreshold (transformation_epsilon_);
  if (transformation_rotation_epsilon_ > 0)
    convergence_criteria_->setRotationThreshold (transformation_rotation_epsilon_);
  else
    convergence_criteria_->setRotationThreshold (1.0 - transformation_epsilon_);

  // Repeat until convergence
  do
  {
    snd_converged = converged_;
    // Get blob data if needed
    PCLPointCloud2::Ptr input_transformed_blob;
    if (need_source_blob_)
    {
      input_transformed_blob.reset (new PCLPointCloud2);
      toPCLPointCloud2 (*input_transformed, *input_transformed_blob);
    }
    // Save the previously estimated transformation
    previous_transformation_ = transformation_;

    // Set the source each iteration, to ensure the dirty flag is updated
    correspondence_estimation_->setInputSource (input_transformed);
    if (correspondence_estimation_->requiresSourceNormals ())
      correspondence_estimation_->setSourceNormals (input_transformed_blob);
    // Estimate correspondences
    if (use_reciprocal_correspondence_)
      correspondence_estimation_->determineReciprocalCorrespondences (*correspondences_, corr_dist_threshold_);
    else
      correspondence_estimation_->determineCorrespondences (*correspondences_, corr_dist_threshold_);

    //if (correspondence_rejectors_.empty ())
    CorrespondencesPtr temp_correspondences (new Correspondences (*correspondences_));
    for (size_t i = 0; i < correspondence_rejectors_.size (); ++i)
    {
      registration::CorrespondenceRejector::Ptr& rej = correspondence_rejectors_[i];
      PCL_DEBUG ("Applying a correspondence rejector method: %s.\n", rej->getClassName ().c_str ());
      if (rej->requiresSourcePoints ())
        rej->setSourcePoints (input_transformed_blob);
      if (rej->requiresSourceNormals () && source_has_normals_)
        rej->setSourceNormals (input_transformed_blob);
      rej->setInputCorrespondences (temp_correspondences);
      rej->getCorrespondences (*correspondences_);
      // Modify input for the next iteration
      if (i < correspondence_rejectors_.size () - 1)
        *temp_correspondences = *correspondences_;
    }

    size_t cnt = correspondences_->size ();
    // Check whether we have enough correspondences
    if (static_cast<int> (cnt) < min_number_correspondences_)
    {
      PCL_ERROR ("[pcl::%s::computeTransformation] Not enough correspondences found. Relax your threshold parameters.\n", getClassName ().c_str ());
      convergence_criteria_->setConvergenceState(pcl::registration::DefaultConvergenceCriteria<Scalar>::CONVERGENCE_CRITERIA_NO_CORRESPONDENCES);
      converged_ = false;
      break;
    }

    // Estimate the transform
    transformation_estimation_->estimateRigidTransformation (*input_transformed, *target_, *correspondences_, transformation_);

    if (nr_iterations_++)
    {
      Vector6 g_k = Matrix42Vector6(transformation_ * final_transformation_);
      err = getMSE(transformation_ * final_transformation_);

      if ((err - prev_err)/prev_err > error_overflow_threshold_)
      {
        u_next = u_k = g.rightCols(1);
        prev_err = std::numeric_limits<double>::max();
        u = u.rightCols(2);
        g = g.rightCols(1);
        f = f.rightCols(1);
      }
      else
      {
        prev_err = err;

        g.conservativeResize(g.rows(),g.cols()+1);
        g.col(g.cols()-1) = g_k;

        Vector6 f_k = g_k - u_k;
        f.conservativeResize(f.rows(),f.cols()+1);
        f.col(f.cols()-1) = f_k;

        u_next = get_next_u(u, g, f);
        u.conservativeResize(u.rows(),u.cols()+1);
        u.col(u.cols()-1) = u_next;

        u_k = u_next;
      }
    }
    else
    {
      Vector6 u0 = Matrix42Vector6(guess);
      u.conservativeResize(u.rows(),u.cols()+1);
      u.col(0)=u0;

      Vector6 u1 = Matrix42Vector6(transformation_ * final_transformation_);
      g.conservativeResize(g.rows(),g.cols()+1);
      g.col(0)=u1;

      u.conservativeResize(u.rows(),u.cols()+1);

      u.col(1)=u1;

      f.conservativeResize(f.rows(),f.cols()+1);
      f.col(0)=u1 - u0;

      u_next = u1;
      u_k = u1;
      prev_err = getMSE(transformation_ * final_transformation_);
    }

    transformation_ = Vector62Matrix4(u_next)*(final_transformation_.inverse());
    final_transformation_ = Vector62Matrix4(u_next);

    // Tranform the data
    transformCloud (*input_, *input_transformed, final_transformation_);


    // Update the vizualization of icp convergence
    //if (update_visualizer_ != 0)
    //  update_visualizer_(output, source_indices_good, *target_, target_indices_good );

    converged_ = static_cast<bool> ((*convergence_criteria_));
  }
  while (!converged_ || !(snd_converged || nr_iterations_ <= small_step_threshold_ || nr_iterations_ == convergence_criteria_->getMaximumIterations()));

  // Transform the input cloud using the final transformation
  PCL_DEBUG ("Transformation is:\n\t%5f\t%5f\t%5f\t%5f\n\t%5f\t%5f\t%5f\t%5f\n\t%5f\t%5f\t%5f\t%5f\n\t%5f\t%5f\t%5f\t%5f\n", 
      final_transformation_ (0, 0), final_transformation_ (0, 1), final_transformation_ (0, 2), final_transformation_ (0, 3),
      final_transformation_ (1, 0), final_transformation_ (1, 1), final_transformation_ (1, 2), final_transformation_ (1, 3),
      final_transformation_ (2, 0), final_transformation_ (2, 1), final_transformation_ (2, 2), final_transformation_ (2, 3),
      final_transformation_ (3, 0), final_transformation_ (3, 1), final_transformation_ (3, 2), final_transformation_ (3, 3));

  // Copy all the values
  output = *input_;
  // Transform the XYZ + normals
  transformCloud (*input_, output, final_transformation_);
}

template <typename PointSource, typename PointTarget, typename Scalar> void
pcl::AndersonIterativeClosestPointWithNormals<PointSource, PointTarget, Scalar>::transformCloud(
  const PointCloudSource &input,
  PointCloudSource &output,
  const Matrix4 &transform)
{
  pcl::transformPointCloudWithNormals(input, output, transform);
}

#endif //PCL_REGISTRATION_IMPL_ICP_ANDERSON_HPP_
